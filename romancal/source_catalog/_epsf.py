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

import matplotlib.pyplot as plt
import numpy as np
import roman_datamodels as rdm
from astropy.nddata import NDData, StdDevUncertainty
from astropy.table import Table
from astropy.time import Time
from astropy.visualization import simple_norm
from photutils.psf import EPSFBuilder, extract_stars
from photutils.psf.flags import PSF_FLAGS
from roman_datamodels.datamodels import EpsfRefModel

from romancal.resample._utils import compute_var_sky
from romancal.source_catalog._calibrate import (
    make_image_mask,
    pixel_scale_arcsec,
    subtract_background,
)
from romancal.source_catalog._utils import copy_model_arrays
from romancal.source_catalog.psf import _get_jitter_params
from romancal.source_catalog.source_catalog_step import (
    SourceCatalogStep,
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


def build_epsf(stars, **kwargs_epsfbuilder):
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
        "--stamp-size",
        type=int,
        default=21,
        help="Size of the extracted stars cutouts. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--reject-core-dq",
        type=int,
        default=int(_DQ.DO_NOT_USE),
        help="""Reject stars with these DQ flags in the central 3x3 pixels.
                Default: %(default)s""",
    )
    pick_star_parser.add_argument(
        "--allow-saturated",
        action="store_true",
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
        "--max-kron-contamination",
        type=float,
        default=0.5,
        help="""Maximum kron_nomask_flux / kron_flux - 1, the fraction of
                neighbor flux in the Kron aperture. Default: %(default)s""",
    )
    pick_star_parser.add_argument(
        "--allow-extended",
        action="store_true",
        help="Keep sources flagged is_extended in the catalog.",
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
        help="""The desired accuracy for the centers of stars. Default: the
                EPSFBuilder default.""",
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

    # Configure logging for all of romancal, as the other scripts do.
    romancal_logger = logging.getLogger("romancal")
    logger_handler = logging.StreamHandler()
    romancal_logger.addHandler(logger_handler)
    logger_format_debug = logging.Formatter(
        "%(levelname)s:%(filename)s::%(funcName)s: %(message)s"
    )
    level = LOGLEVELS[min(len(LOGLEVELS) - 1, args.verbose)]
    romancal_logger.setLevel(level)
    if level <= logging.DEBUG:
        logger_handler.setFormatter(logger_format_debug)
    logger.info("epsf_from_model called with args %s", args)

    kwargs_pick_stars = groups["pick_star"]
    kwargs_pick_stars["reject_saturated"] = not kwargs_pick_stars.pop("allow_saturated")
    kwargs_pick_stars["reject_extended"] = not kwargs_pick_stars.pop("allow_extended")

    # Determine ePSF
    with rdm.open(args.model) as model:
        catalog = Table.read(args.catalog)
        epsf_model, epsf, stars = epsf_from_model(
            model,
            catalog,
            kwargs_pick_stars=kwargs_pick_stars,
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
        The jitter defaults to the image's guide star jitter, since the
        ePSF is measured from stars that include it.

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
    stars = pick_stars(model, catalog, **kwargs_pick_stars)

    if kwargs_build_epsf is None:
        kwargs_build_epsf = dict()
    epsf = build_epsf(stars, **kwargs_build_epsf)

    if meta is None:
        meta = dict()
    jitter = _get_jitter_params(getattr(model.meta, "guide_star", dict()))
    meta = jitter | meta
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
    stamp_size=21,
    reject_core_dq=_DQ.DO_NOT_USE,
    reject_saturated=True,
    max_masked_fraction=0.2,
    max_psf_cfit=0.2,
    min_snr=50.0,
    max_kron_contamination=0.5,
    reject_extended=True,
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

    reject_core_dq : int or None
        Stars with any of these DQ flags in the central 3x3 pixels are
        rejected.

    reject_saturated : bool
        Reject stars with any SATURATED pixel in their stamp. Saturated
        pixels are measured from early resultants only and so see
        different jitter than unsaturated pixels.

    max_masked_fraction : float or None
        Reject stars with more than this fraction of their stamp masked.

    max_psf_cfit : float or None
        Maximum catalog ``psf_cfit``, the PSF fit residual in the central
        pixel divided by the fitted flux. Stars are below about 0.1;
        single hot pixels and cosmic rays are typically above 1.

    min_snr : float or None
        Minimum ``psf_flux / psf_flux_err``.

    max_kron_contamination : float or None
        Maximum ``kron_nomask_flux / kron_flux - 1``, the flux from
        neighboring sources in the Kron aperture relative to the source's
        own Kron flux. Fainter neighbors in the wings of larger stamps
        are removed by the sigma clipping in EPSFBuilder.

    reject_extended : bool
        Reject sources flagged ``is_extended`` in the catalog.

    reject_psf_flags : int or None
        Reject objects with any of these ``psf_flags`` bits set.

    bkg_boxsize : float or None
        Background mesh box size in arcsec, as in the source catalog
        step. If None, use the source catalog step default.

    Returns
    -------
    stars : EPSFStars
        The stars to use for ePSF determination.
    """
    logger.info("Picking stars from catalog of %s objects.", len(catalog))

    # Mask and subtract the background as the source catalog step does,
    # also masking pixels without a usable noise estimate. Data in the
    # masked pixels are never used. The noise is the sky noise used by
    # resample, which, unlike err, does not depend on the source flux;
    # weights that drop with flux would make the star fits depend on
    # the flux.
    dq = np.asarray(model.dq)
    with np.errstate(invalid="ignore"):
        err = np.sqrt(compute_var_sky(model))
    pixmask = make_image_mask(model) | ~np.isfinite(err) | (err <= 0)
    if bkg_boxsize is None:
        bkg_boxsize = SourceCatalogStep().bkg_boxsize
    bkg_boxsize_px = max(round(bkg_boxsize / pixel_scale_arcsec(model)), 1)
    model = copy_model_arrays(model)
    subtract_background(model, pixmask, bkg_boxsize_px)
    data = np.where(pixmask, 0.0, model.data)
    err = np.where(pixmask, 1.0, err)

    # Catalog level cuts. Sources without a PSF fit position are dropped.
    x = np.asarray(catalog["x_psf"], dtype=float)
    y = np.asarray(catalog["y_psf"], dtype=float)
    hsize = (stamp_size - 1) / 2
    ny, nx = data.shape
    cuts = {
        "in image": (
            (x > hsize) & (x < nx - 1 - hsize) & (y > hsize) & (y < ny - 1 - hsize)
        ),
    }
    if reject_psf_flags is not None:
        cuts["psf fit"] = (catalog["psf_flags"] & reject_psf_flags) == 0
    if reject_extended:
        cuts["point source"] = ~catalog["is_extended"]
    if max_psf_cfit is not None:
        cuts["psf cfit"] = catalog["psf_cfit"] <= max_psf_cfit
    if min_snr is not None:
        cuts["snr"] = catalog["psf_flux"] / catalog["psf_flux_err"] >= min_snr
    if max_kron_contamination is not None:
        kron_ratio = catalog["kron_nomask_flux"] / catalog["kron_flux"]
        cuts["kron contamination"] = kron_ratio - 1 <= max_kron_contamination
    keep = np.ones(len(catalog), dtype=bool)
    for name, cut in cuts.items():
        keep &= np.asarray(cut)
        logger.debug("After %s cut: %s", name, keep.sum())

    # Pixel level cuts
    sat_flag = _DQ.SATURATED
    for idx in np.flatnonzero(keep):
        ix, iy = round(x[idx]), round(y[idx])
        h = int(hsize)
        stamp = (slice(iy - h, iy + h + 1), slice(ix - h, ix + h + 1))
        core = (slice(iy - 1, iy + 2), slice(ix - 1, ix + 2))
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

    # Extract stars
    stars_tbl = Table()
    stars_tbl["x"] = x[keep]
    stars_tbl["y"] = y[keep]
    stars_tbl["id"] = np.asarray(catalog["label"])[keep]
    nddata = NDData(data, mask=pixmask, uncertainty=StdDevUncertainty(err))
    stars = extract_stars(nddata, stars_tbl, size=stamp_size)

    # That's all folks
    logger.info("Found %s valid stars", len(stars))
    return stars


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
