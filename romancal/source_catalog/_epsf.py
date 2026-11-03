"""Effective Point Spread Function utilities.

Utilities to create ePSF models to be used where PSFs are needed in source detection
and flux/position analysis.

`epsf_from_model` encapsulates the general order of operations which is:
- pick_stars
- build_epsf
- make_epsf_reference

`cli` is the command-line interface entrypoint.
"""

import argparse
import logging
import math
import random
import warnings

import matplotlib.pyplot as plt
import numpy as np
import roman_datamodels as rdm
from astropy.nddata import NDData, StdDevUncertainty
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time
from astropy.visualization import simple_norm
from photutils.psf import EPSFBuilder, extract_stars
from photutils.psf.flags import PSF_FLAGS
from roman_datamodels.datamodels import EpsfRefModel
from scipy.spatial import cKDTree

from romancal.resample._utils import compute_var_sky
from romancal.source_catalog._background import RomanBackground
from romancal.source_catalog.source_catalog_step import (
    SourceCatalogStep,
    _pixel_scale,
)

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
DEBUG_FULL = logging.DEBUG - 1
LOGLEVELS = [logging.INFO, logging.DEBUG, DEBUG_FULL]
_DQ = rdm.dqflags.pixel

# psf_flags bits that make a PSF fit unusable for picking ePSF stars.
# N_PIXELS_FIT_PARTIAL is left out; masked pixels are handled by pick_stars.
REJECT_PSF_FLAGS = (
    PSF_FLAGS.POSITION_OUTSIDE_IMAGE
    | PSF_FLAGS.NON_POSITIVE_FLUX
    | PSF_FLAGS.NO_CONVERGENCE
    | PSF_FLAGS.NO_COVARIANCE
    | PSF_FLAGS.PARAMETER_NEAR_BOUND
    | PSF_FLAGS.NO_OVERLAP
    | PSF_FLAGS.FULLY_MASKED
    | PSF_FLAGS.TOO_FEW_PIXELS
    | PSF_FLAGS.NON_FINITE_POSITION
    | PSF_FLAGS.NON_FINITE_FLUX
    | PSF_FLAGS.POSITION_OUTSIDE_FIT_REGION
)


def build_epsf(stars, max_psferr=None, **kwargs_epsfbuilder):
    """Create the ePSF from a list of stars.

    Parameters
    ----------
    stars : EPSFStars
        The stars to base the ePSF on.

    kwargs_epsfbuilder : dict
        Arguments to pass to EPSFBuilder. If the following arguments
        are left unspecified, these defaults will be used:
        - oversample=4
        - maxiters=30
        - center_accuracy: Value of argument `max_psferr` if max_psferr is not None.
        - progress_bar=false

    Returns
    -------
    epsf_result : EPFSBuildResults
        The created ePSF.
    """
    logger.info("Building ePSF from catalog of %s stars", len(stars))
    epsfbuilder_final_kwargs = {
        "oversampling": 4,
        "maxiters": 30,
        "progress_bar": False,
    }
    if kwargs_epsfbuilder is None:
        kwargs_epsfbuilder = dict()
    epsfbuilder_final_kwargs.update(kwargs_epsfbuilder)
    if epsfbuilder_final_kwargs.get("center_accuracy") is None:
        try:
            del epsfbuilder_final_kwargs["center_accuracy"]
        except KeyError:
            pass
    if max_psferr is not None:
        epsfbuilder_final_kwargs["center_accuracy"] = epsfbuilder_final_kwargs.get(
            "center_accuracy", max_psferr
        )

    epsfbuilder = EPSFBuilder(**epsfbuilder_final_kwargs)
    epsf_result = epsfbuilder(stars)

    # That's all folks.
    logger.info("Calculated ePSF parameters: %s", epsf_result)
    return epsf_result


def cli(argv=None):
    """Command-line interface entry point.

    Given an Roman L2 ImageModel and associated source catalog, create a ePSF reference file.
    """
    parser = argparse.ArgumentParser(
        description="Create an ePSF reference file from a Roman L2 image and associated catalog"
    )

    # General arguments
    parser.add_argument(
        "model", type=str, help="Path to the Roman L2 calibration image."
    )
    parser.add_argument(
        "catalog",
        type=str,
        help="Path to the source catalog associated with the L2 image.",
    )
    parser.add_argument(
        "output_epsf_reffile",
        type=str,
        help="Path to write the newly created ePSF reference file.",
    )
    parser.add_argument("--plot-epsf", action="store_true", help="Plot the ePSF.")
    parser.add_argument(
        "--plot-stamps",
        type=int,
        default=None,
        help="Plot the specified number of randomly selected stamps.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity. Specifying multiple times adds more output.",
    )

    # `pick_star` arguments
    pick_star_parser = parser.add_argument_group("pick_star", "Star picking options.")
    pick_star_parser.add_argument(
        "--stamp-size", type=int, default=9, help="Size of the extracted stars cutouts."
    )
    pick_star_parser.add_argument(
        "--mask-dq",
        type=int,
        default=int(_DQ.DO_NOT_USE),
        help="DQ flags of pixels to mask in the stamps. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--reject-core-dq",
        type=int,
        default=int(_DQ.DO_NOT_USE),
        help="Reject stars with these DQ flags in their core. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--core-size",
        type=int,
        default=3,
        help="Size of the core used by --reject-core-dq. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--allow-saturated",
        action="store_false",
        dest="reject_saturated",
        help="Keep stars with saturated pixels in their stamps.",
    )
    pick_star_parser.add_argument(
        "--max-psf-cfit",
        type=float,
        default=0.2,
        help="""Maximum catalog psf_cfit; rejects hot pixels and cosmic rays.
                Default: %(default)s""",
    )
    pick_star_parser.add_argument(
        "--max-masked-fraction",
        type=float,
        default=0.2,
        help="Maximum masked fraction of a stamp. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--min-snr",
        type=float,
        default=50.0,
        help="Minimum psf_flux / psf_flux_err. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--max-contamination",
        type=float,
        default=0.05,
        help="""Maximum summed neighbor flux within --contamination-radius,
                relative to the star flux. Default: %(default)s""",
    )
    pick_star_parser.add_argument(
        "--contamination-radius",
        type=float,
        default=6.5,
        help="Radius in pixels for --max-contamination. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--max-psferr",
        type=float,
        default=None,
        help="""Maximum  sqrt(x_psf_err**2 + y_psf_err**2) allowed.
                                  If None, sigma-clipped stats to remove anything larger than
                                  two times the deviation of the psf errors for all objects.
                                  Default: %(default)s""",
    )
    pick_star_parser.add_argument(
        "--allow-extended",
        action="store_false",
        dest="mask_extended",
        help="""Allow extended sources.""",
    )
    pick_star_parser.add_argument(
        "--bkg-boxsize",
        type=float,
        default=None,
        help="""Background mesh box size in arcsec. Default: the source
                catalog step default.""",
    )
    pick_star_parser.add_argument(
        "--reject-psf-flags",
        type=int,
        default=int(REJECT_PSF_FLAGS),
        help="Reject objects with these psf_flags bits set. Default: %(default)s",
    )

    # `build_epsf` arguments
    build_epsf_parser = parser.add_argument_group(
        "build_epsf", "photutils.EPSFBuilder options."
    )
    build_epsf_parser.add_argument(
        "--oversampling",
        type=int,
        default=4,
        help="""The integer oversampling factor(s) of the output ePSF relative to the input stars.
                                   Default: %(default)s""",
    )
    build_epsf_parser.add_argument(
        "--maxiters",
        type=int,
        default=30,
        help="""The maximum number of ePSF building iterations to perform.
                                   Default: %(default)s""",
    )
    build_epsf_parser.add_argument(
        "--center-accuracy",
        type=float,
        default=None,
        help="""The desired accuracy for the centers of stars. If not specified, this will be
                                   based on the value of `max_psferr`
                                   Default: %(default)s""",
    )

    args = parser.parse_args(argv)

    # Separate out the group arguments.
    groups = {}
    for group in parser._action_groups:
        if group.title in ["positional arguments", "options"]:
            continue
        group_args = {
            action.dest: getattr(args, action.dest)
            for action in group._group_actions
            if action.dest != "help"
        }
        groups[group.title] = group_args

    # Configure logging
    logger_handler = logging.StreamHandler()
    logger.addHandler(logger_handler)
    logger_format_debug = logging.Formatter(
        "%(levelname)s:%(filename)s::%(funcName)s: %(message)s"
    )
    level = LOGLEVELS[min(len(LOGLEVELS) - 1, args.verbose)]
    logger.setLevel(level)
    if level > logging.DEBUG:
        warnings.simplefilter("ignore", UserWarning)
    else:
        logger_handler.setFormatter(logger_format_debug)
    logger.info("epsf_from_model called with args %s", args)

    # Determine ePSF
    with rdm.open(args.model) as model:
        catalog = Table.read(args.catalog)
        epsf_model, epsf, stars = epsf_from_model(
            model,
            catalog,
            kwargs_pick_stars=groups["pick_star"],
            kwargs_build_epsf=groups["build_epsf"],
        )

    if args.plot_epsf:
        plot_epsfstar(epsf.epsf.data, finalize=True)
    if args.plot_stamps:
        to_plot = random.choices(stars, k=args.plot_stamps)
        plot_epsfstars(to_plot)

    # That's all folks.
    logger.info("Saving ePSF to %s", args.output_epsf_reffile)
    epsf_model.save(args.output_epsf_reffile)


def epsf_from_model(
    model, catalog, kwargs_pick_stars=None, kwargs_build_epsf=None, meta=None
):
    """Given a Roman L2 ImageModel and associated source catalog, construct the ePSF.

    Parameters
    ----------
    model : roman_datamodels.datamodels.ImageModel
        L2 model

    catalog : astropy.table.Table
        catalog based on model

    kwargs_pick_stars : dict or None
        Keyword arguments to pass to the `pick_stars` function.

    kwargs_build_epsf : dict or None
        Keyword arguments to pass to the `build_epsf` function.

    meta : dict or None
        The meta to use to populate the ePSF reference model meta with.

    Returns
    -------
    epsf_model, epsf, stars : EpsfRefModel, ImagePSF, EPSFStars
        The EPSF reference model, the ePSF, and the stars used.
    """
    # Retrieve some basic info from the input model.
    detector = model.meta.instrument.detector
    optical_element = model.meta.instrument.optical_element

    if kwargs_pick_stars is None:
        kwargs_pick_stars = dict()
    stars, max_psferr = pick_stars(model, catalog, **kwargs_pick_stars)

    if kwargs_build_epsf is None:
        kwargs_build_epsf = dict()
    epsf = build_epsf(stars, max_psferr, **kwargs_build_epsf)

    if meta is None:
        meta = dict()
    epsf_model = make_epsf_reference(epsf, detector, optical_element, meta=meta)

    # That's all folks
    return epsf_model, epsf, stars


def make_epsf_reference(epsf, detector, optical_element, meta=None):
    """Create an EPSF reference model from an ePSF.

    Parameters
    ----------
    epsf : EPSFBuildResult
        The epsf to store.

    detector : str
        The detector represented. Should be of form "WFIXX"
        where "XX" is zero-padded detector number

    optical_element : str
        The optical_element represented.

    meta : dict
        Reference meta values to use. If unspecified, minimal fake values will be used.

    Returns
    -------
    epsf_model : EpsfRefModel
        The EPSF Reference model.
    """
    logger.info("Creating ePSF reference.")
    default_meta = {
        "author": "make_epsf_reference",
        "description": str(epsf),
        "jitter_major": 0.0,
        "jitter_minor": 0.0,
        "jitter_position_angle": 0.0,
        "instrument": {
            "name": "WFI",
            "detector": detector.upper(),
            "optical_element": optical_element.upper(),
        },
        "origin": "STSCI",
        "oversample": epsf.epsf.oversampling[0],
        "pedigree": f"INFLIGHT {Time.now()}",
        "pixel_x": [2047.5],
        "pixel_y": [2047.5],
        "telescope": "ROMAN",
        "useafter": Time("2026-09-30", format="isot"),
    }
    if meta is None:
        meta = dict()
    final_meta = default_meta | meta
    defaults = {"meta": final_meta}

    epsf_data = epsf.epsf.data
    psf = np.zeros((3, 3, 1, epsf_data.shape[0], epsf_data.shape[1]), dtype=np.float32)
    psf[0, 1, 0, :, :] = epsf_data
    defaults["psf"] = psf
    model = EpsfRefModel.create_minimal(defaults=defaults)

    # Update meta that are not defined in the model.
    # These need be done after model creation or they are ignored.
    model.meta["oversample"] = meta.get("oversample", epsf.epsf.oversampling[0])

    return model


def pick_stars(
    model,
    catalog,
    stamp_size=9,
    mask_dq=_DQ.DO_NOT_USE,
    reject_core_dq=_DQ.DO_NOT_USE,
    core_size=3,
    reject_saturated=True,
    max_masked_fraction=0.2,
    max_psf_cfit=0.2,
    min_snr=50.0,
    max_contamination=0.05,
    contamination_radius=6.5,
    max_psferr=None,
    mask_extended=True,
    reject_psf_flags=REJECT_PSF_FLAGS,
    bkg_boxsize=None,
):
    """Pick stars out of an image based on various criteria.

    Parameters
    ----------
    model : roman_datamodels.datamodels.ImageModel
        L2 model

    catalog : astropy.table.Table
        catalog based on model

    stamp_size : int
        Size of the extracted stars cutouts.

    mask_dq : int or None
        Pixels with any of these DQ flags are masked (not used) in the
        stamps. Non-finite pixels are always masked.

    reject_core_dq : int or None
        Stars with any of these DQ flags in the central
        ``core_size`` x ``core_size`` pixels are rejected.

    core_size : int
        Size of the core region used by ``reject_core_dq``.

    reject_saturated : bool
        Reject stars with any SATURATED pixel in their stamp. Saturated
        pixels are measured from early resultants only and so see
        different jitter than unsaturated pixels.

    max_masked_fraction : float or None
        Reject stars with more than this fraction of their stamp masked.

    max_psf_cfit : float or None
        Maximum catalog ``psf_cfit``, the PSF fit residual in the central
        pixel divided by the fitted flux. Stars are below about 0.1;
        single hot pixels and cosmic rays are typically above 1. Not
        applied if the catalog has no ``psf_cfit`` column.

    min_snr : float or None
        Minimum ``psf_flux / psf_flux_err``.

    max_contamination : float or None
        Maximum summed flux of catalog neighbors within
        ``contamination_radius``, relative to the star's flux. Fluxes are
        taken from ``segment_flux``.

    contamination_radius : float
        Radius in pixels used for ``max_contamination``. It does not
        grow with ``stamp_size``; fainter neighbors in the wings of
        larger stamps are removed by the sigma clipping in EPSFBuilder.

    max_psferr : float or None
        Maximum  sqrt(x_psf_err**2 + y_psf_err**2) allowed.
        If None, sigma-clipped stats to remove anything larger than
        the median plus two times the deviation of the psf errors.

    mask_extended : bool
        Remove extended sources.

    reject_psf_flags : int or None
        Reject objects with any of these ``psf_flags`` bits set, or with
        non-finite PSF position errors.

    bkg_boxsize : float or None
        Background mesh box size in arcsec, as in the source catalog
        step. If None, use the source catalog step default.

    Returns
    -------
    stars : EPSFStars
        The stars to use for ePSF determination.

    max_psferr : float
        The maximum psf error used.
    """
    logger.info("Picking stars from catalog of %s objects.", len(catalog))

    # Background; data in the masked pixels are never used.
    dq = np.asarray(model.dq)
    err = background_noise(model)
    data = np.asarray(model.data, dtype=float)
    pixmask = ~np.isfinite(data) | ~np.isfinite(model.err)
    pixmask |= ~np.isfinite(err) | (err <= 0)
    if mask_dq is not None:
        pixmask |= (dq & mask_dq) > 0
    if bkg_boxsize is None:
        bkg_boxsize = SourceCatalogStep().bkg_boxsize
    bkg_boxsize_px = max(round(bkg_boxsize / _pixel_scale(model)), 1)
    background = RomanBackground(data, bkg_boxsize_px, mask=pixmask).background
    data = np.where(pixmask, 0.0, data - background)
    err = np.where(pixmask, 1.0, err)

    # Catalog level cuts
    x = _filled(catalog["x_psf"])
    y = _filled(catalog["y_psf"])
    use_centroid = ~np.isfinite(x) | ~np.isfinite(y)
    x[use_centroid] = _filled(catalog["x_centroid"])[use_centroid]
    y[use_centroid] = _filled(catalog["y_centroid"])[use_centroid]
    psf_err = np.hypot(_filled(catalog["x_psf_err"]), _filled(catalog["y_psf_err"]))
    hsize = (stamp_size - 1) / 2
    ny, nx = data.shape
    cuts = {
        "in image": (
            (x > hsize) & (x < nx - 1 - hsize) & (y > hsize) & (y < ny - 1 - hsize)
        ),
    }
    if reject_psf_flags is not None:
        psf_flags = np.ma.filled(np.ma.asarray(catalog["psf_flags"]), reject_psf_flags)
        cuts["psf fit"] = ((psf_flags & reject_psf_flags) == 0) & np.isfinite(psf_err)
    if mask_extended:
        cuts["point source"] = ~_filled(catalog["is_extended"], True).astype(bool)
    if max_psf_cfit is not None:
        if "psf_cfit" in catalog.colnames:
            cuts["psf cfit"] = ~(_filled(catalog["psf_cfit"]) > max_psf_cfit)
        else:
            logger.info("Catalog has no psf_cfit; not cutting on it.")
    if min_snr is not None:
        snr = _filled(catalog["psf_flux"]) / _filled(catalog["psf_flux_err"])
        cuts["snr"] = snr >= min_snr
    if max_contamination is not None:
        cuts["contamination"] = (
            neighbor_contamination(catalog, x, y, contamination_radius)
            <= max_contamination
        )
    keep = np.ones(len(catalog), dtype=bool)
    for name, cut in cuts.items():
        keep &= cut
        logger.debug("After %s cut: %s", name, keep.sum())

    # Pixel level cuts
    sat_flag = _DQ.SATURATED
    hcore = core_size // 2
    for idx in np.flatnonzero(keep):
        ix, iy = round(x[idx]), round(y[idx])
        h = int(hsize)
        stamp = (slice(iy - h, iy + h + 1), slice(ix - h, ix + h + 1))
        core = (slice(iy - hcore, iy + hcore + 1), slice(ix - hcore, ix + hcore + 1))
        if pixmask[core].all():
            keep[idx] = False
        elif reject_core_dq is not None and (dq[core] & reject_core_dq).any():
            keep[idx] = False
        elif reject_saturated and (dq[stamp] & sat_flag).any():
            keep[idx] = False
        elif max_masked_fraction is not None and (
            pixmask[stamp].mean() > max_masked_fraction
        ):
            keep[idx] = False
    logger.debug("After pixel cuts: %s", keep.sum())

    # With all gross undesirable detections removed,
    # make a final sweep using the psf error.
    if max_psferr is None:
        _, median, stddev = sigma_clipped_stats(psf_err[keep])
        max_psferr = median + (2 * stddev)
    keep &= psf_err <= max_psferr

    # Extract stars
    stars_tbl = Table()
    stars_tbl["x"] = x[keep]
    stars_tbl["y"] = y[keep]
    stars_tbl["id"] = np.asarray(catalog["label"])[keep]
    nddata = NDData(data, mask=pixmask, uncertainty=StdDevUncertainty(err))
    stars = extract_stars(nddata, stars_tbl, size=stamp_size)

    # That's all folks
    logger.info(
        "Found %s valid stars with maximum PSF error of %s", len(stars), max_psferr
    )
    return stars, max_psferr


def background_noise(model):
    """Per-pixel noise that does not depend on the source flux.

    Uncertainties that include source Poisson noise give bright pixels
    low weight, which makes the star fits depend on the flux. Instead
    use the sky variance used by resample: read noise plus the Poisson
    variance scaled to the median sky level.

    Parameters
    ----------
    model : roman_datamodels.datamodels.ImageModel
        L2 model

    Returns
    -------
    sigma : np.ndarray
        Standard deviation per pixel; NaN where the sky variance is not
        positive.
    """
    var_sky = np.asarray(compute_var_sky(model), dtype=float)
    with np.errstate(invalid="ignore"):
        return np.where(var_sky > 0, np.sqrt(var_sky), np.nan)


def neighbor_contamination(catalog, x, y, radius, flux_col="segment_flux"):
    """Summed flux of catalog neighbors relative to each source's flux.

    Parameters
    ----------
    catalog : astropy.table.Table
        Source catalog.

    x, y : np.ndarray
        Source positions in pixels.

    radius : float
        Neighbors within this distance in pixels are counted.

    flux_col : str
        Catalog column with the flux of each source.

    Returns
    -------
    contamination : np.ndarray
        Sum of neighbor fluxes divided by the source flux. Sources with
        non-positive or non-finite flux are assigned inf.
    """
    flux = _filled(catalog[flux_col])
    good = np.isfinite(x) & np.isfinite(y)
    pos = np.column_stack([x, y])
    tree = cKDTree(pos[good])
    gflux = np.clip(np.nan_to_num(flux[good]), 0, None)
    contamination = np.full(len(catalog), np.inf)
    gidx = np.flatnonzero(good)
    neighbor_lists = tree.query_ball_point(pos[good], radius)
    for k, (i, neighbors) in enumerate(zip(gidx, neighbor_lists, strict=True)):
        if flux[i] > 0:
            contamination[i] = (gflux[neighbors].sum() - gflux[k]) / flux[i]
    return contamination


def _filled(column, fill_value=np.nan):
    """Return a catalog column as a plain array, filling masked values."""
    values = np.ma.filled(np.ma.asarray(column), fill_value)
    return np.asarray(getattr(values, "value", values), dtype=float)


def plot_epsfstar(star, ax=None, title=None, finalize=False, **kwargs_simple_norm):
    """Plot EPSFStar.

    Parameters
    ----------
    star : EPSFStar or other array-like object
        The 2D data to plot

    ax : matplotlib.Axes or None
        The axis to plot on. If None, use the current axis.

    finalize : bool
        Close out the plot.

    kwargs_simple_norm : dict
        Arguments for the `astropy.visualization.simple_norm` function. If None, the "log" stretch with "percent=99.0" is used.
    """
    if not len(kwargs_simple_norm):
        kwargs_simple_norm = {"stretch": "log", "percent": 99.0}
    norm = simple_norm(star, **kwargs_simple_norm)
    if ax is None:
        ax = plt.gca()
    if title is not None:
        ax.set_title(title)
    ax.imshow(star, norm=norm, origin="lower", cmap="viridis")
    if finalize:
        plt.show()


def plot_epsfstars(stars, ncols=5, title_func=None, **kwargs_simple_norm):
    """Plot EPSFStars

    Parameters
    ----------
    stars : EPSFStars
        stars to plot

    ncols : int
        Number of columns

    title_func : func
        Function to generate a string to use as the title for each star.
        Expect arguments `(axis, idx, star)`.

    kwargs_simple_norm : dict
        Arguments for the `simple_norm` function. If None, the "log" stretch with "percent=99.0" is used.
    """
    if title_func is None:
        title_func = lambda axis, idx, star: f"Star {idx}"

    nrows = math.ceil(len(stars) / ncols)
    fig, ax = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4 * ncols, 4 * nrows),
        squeeze=False,
        layout="compressed",
    )
    ax = ax.ravel()

    for idx, (axis, star) in enumerate(zip(ax, stars, strict=False)):
        plot_epsfstar(
            star, ax=axis, title=title_func(axis, idx, star), **kwargs_simple_norm
        )

    plt.show()
