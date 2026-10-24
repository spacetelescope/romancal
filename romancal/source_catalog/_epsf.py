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

import numpy as np
import roman_datamodels as rdm
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
from photutils.psf import EPSFBuilder, EPSFStars, extract_stars
from roman_datamodels.datamodels import EpsfRefModel

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
DEBUG_FULL = logging.DEBUG - 1
LOGLEVELS = [logging.INFO, logging.DEBUG, DEBUG_FULL]


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
    parser.add_argument(
        '--plot-epsf', action="store_true",
        help='Plot the ePSF.'
    )
    parser.add_argument(
        '--plot-stamps', type=int, default=None,
        help='Plot the specified number of randomly selected stamps.'
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
        "--allow-nans",
        action="store_false",
        dest="mask_nans",
        help="Allows stars with NaNs",
    )
    pick_star_parser.add_argument(
        "--mask-dq",
        type=int,
        default=(rdm.dqflags.pixel.DO_NOT_USE + rdm.dqflags.pixel.SATURATED),
        help="Pixel data quality flags used to mask questionable pixels from the data. Default: %(default)s",
    )
    pick_star_parser.add_argument(
        "--min-flux",
        type=float,
        default=None,
        help="""Minimum flux to consider.
                                  Flux is determined by the `photoutils.psf.extract_stars` function. Default: %(default)s""",
    )
    pick_star_parser.add_argument(
        "--min-stdev",
        type=float,
        default=0.3,
        help="""Minimum standard deviation with a stamp to be valid.
                                  The stats are calculated using sigma-clipped determination over each object's stamp.
                                  Default: %(default)s""",
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
        "--allow-failed-psffit",
        action="store_false",
        dest="mask_failed_psffit",
        help="""Removed objects that have a bad PSF fit. Indicated by the x/y_psf_err being 0.""",
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
    mask_nans=True,
    mask_dq=(rdm.dqflags.pixel.DO_NOT_USE + rdm.dqflags.pixel.SATURATED),
    min_flux=500.0,
    min_stdev=None,
    max_psferr=None,
    mask_extended=True,
    mask_failed_psffit=True,
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

    mask_nans : bool
        Remove stars that have NaNs.

    mask_dq : int or None
        Pixel data quality flags used to mask questionable pixels from the data.

    min_flux : float or None
        Minimum flux to consider. Flux is determined by the `photoutils.psf.extract_stars` function.

    min_stdev : float or None
        Minimum standard deviation with a stamp to be valid.
        The stats are calculated using sigma-clipped determination over each object's stamp.

    max_psferr : float or None
        Maximum  sqrt(x_psf_err**2 + y_psf_err**2) allowed.
        If None, sigma-clipped stats to remove anything larger than
        two times the deviation of the psf errors for all objects.

    mask_extended : bool
        Remove extended sources.

    mask_failed_psffit : bool
        Removed objects that have a bad PSF fit. Indicated by the x/y_psf_err being 0.

    Returns
    -------
    stars : EPSFStars
        The stars to use for ePSF determination.

    max_psferr : float
        Two * deviation of the individual stars' psf error estimates, determined using sigma-clipped statistics.
    """
    logger.info("Picking stars from catalog of %s objects.", len(catalog))

    # Massage the data
    mean_val, median_val, std_val = sigma_clipped_stats(model.data, sigma=2.0)
    data = model.data - median_val

    # Create simplified table
    hsize = (stamp_size - 1) / 2
    x = catalog["x_centroid"]
    y = catalog["y_centroid"]
    is_extended = catalog["is_extended"]
    label = catalog["label"]
    psf_err = np.sqrt(catalog["x_psf_err"] ** 2 + catalog["y_psf_err"] ** 2)
    mask = (
        (x > hsize)
        & (x < (data.shape[1] - 1 - hsize))
        & (y > hsize)
        & (y < (data.shape[0] - 1 - hsize))
    )
    if mask_extended:
        mask = [m & (not ie) for m, ie in zip(mask, is_extended, strict=True)]
    if mask_failed_psffit:
        mask = [m & (psferr != 0.0) for m, psferr in zip(mask, psf_err, strict=True)]

    stars_tbl = Table()
    stars_tbl["x"] = x[mask]
    stars_tbl["y"] = y[mask]
    stars_tbl["label"] = label[mask]
    stars_tbl["psf_err"] = psf_err[mask]

    # Extract stars
    stars = extract_stars(NDData(data), stars_tbl, size=stamp_size)
    for star, label, psf_err in zip(
        stars, stars_tbl["label"], stars_tbl["psf_err"], strict=True
    ):
        star.label = label
        star.psf_err = psf_err
        star.stats = sigma_clipped_stats(star.data, sigma=2.0)

    # Throw out stars with bad values
    if mask_nans:
        ignore_mask = [np.any(np.isnan(star)) for star in stars]
    else:
        ignore_mask = [False] * len(stars)

    # Throw out stars with bad dq flags
    if mask_dq is not None:

        def dq_mask(dq, star, mask_dq):
            bb = star.bbox
            dq = dq[bb.iymin : bb.iymax, bb.ixmin : bb.ixmax]
            mask = (dq & mask_dq).any()
            return mask

        ignore_mask = [
            mask or dq_mask(model.dq, star, mask_dq)
            for mask, star in zip(ignore_mask, stars, strict=True)
        ]

    # Only bright
    if min_flux is not None:
        ignore_mask = [
            mask or star.estimate_flux() < min_flux
            for mask, star in zip(ignore_mask, stars, strict=True)
        ]

    # Only significant detections
    if min_stdev is not None:
        ignore_mask = [
            mask or star.stats[2] < min_stdev
            for mask, star in zip(ignore_mask, stars, strict=True)
        ]

    # With all gross undesirable detections removed,
    # make a final sweep using the psf error.
    if max_psferr is None:
        all_psf_errs = [
            star.psf_err
            for mask, star in zip(ignore_mask, stars, strict=True)
            if not mask
        ]
        _, median, stddev = sigma_clipped_stats(all_psf_errs)
        max_psferr = median + (2 * stddev)
    ignore_mask = [
        mask or star.psf_err > max_psferr
        for mask, star in zip(ignore_mask, stars, strict=True)
    ]

    # Filter the stars based on the ignore mask
    stars = EPSFStars(
        [star for mask, star in zip(ignore_mask, stars, strict=True) if not mask]
    )

    # That's all folks
    logger.info(
        "Found %s valid stars with maximum PSF error of %s", len(stars), max_psferr
    )
    return stars, max_psferr


def plot_epsfstar(star, ax=plt, title=None, finalize=False, **kwargs_simple_norm):
    """Plot EPSFStar.

    Parameters
    ----------
    star : EPSFStar or other array-like object
        The 2D data to plot

    ax : matplotlib.Axes
        The axis to plot on.

    finalize : bool
        Close out the plot.

    kwargs_simple_norm : dict
        Arguments for the `astropy.visualization.simple_norm` function. If None, the "log" stretch with "percent=99.0" is used.
    """
    if not len(kwargs_simple_norm):
        kwargs_simple_norm = {'stretch': 'log', 'percent': 99.0}
    norm = simple_norm(star, **kwargs_simple_norm)
    if title is not None:
        ax.set_title(title)
    ax.imshow(star, norm=norm, origin='lower', cmap='viridis')
    if finalize:
        ax.show()


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
        title_func = lambda axis, idx, star: f'Star {idx}'

    nstars = len(stars)
    nrows = math.ceil(len(stars) / ncols)
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 4 * nstars),
                           squeeze=True, layout='compressed',
                          )
    if ncols > 1:
        ax = ax.ravel()

    for idx, (axis, star) in enumerate(zip(ax, stars)):
        plot_epsfstar(star, ax=axis, title=title_func(axis, idx, star), **kwargs_simple_norm)

    plt.show()
