"""
Utilities for fitting model PSFs to rate images.
"""

import logging
import os

import astropy.units as u
import numpy as np
import webbpsf
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import CCDData, bitmask
from astropy.table import Table
from photutils.background import LocalBackground
from photutils.detection import DAOStarFinder
from photutils.psf import (
    GriddedPSFModel,
    IterativePSFPhotometry,
    PSFPhotometry,
    SourceGrouper,
)
from roman_datamodels.datamodels import ImageModel
from webbpsf import conf, gridded_library, restart_logging

from romancal.lib.dqflags import pixel as roman_dq_flag_map

__all__ = [
    "create_gridded_psf_model",
    "fit_psf_to_image_model",
    "dq_to_boolean_mask",
]

# set loggers to debug level by default:
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Phase C central wavelengths [micron], released by Goddard (Jan 2023):
# https://roman.ipac.caltech.edu/sims/Param_db.html#wfi_filters
filter_central_wavelengths = {
    "F062": 0.620,
    "F087": 0.869,
    "F106": 1.060,
    "F129": 1.293,
    "F146": 1.464,
    "F158": 1.577,
    "F184": 1.842,
    "F213": 2.125,
}

default_finder = DAOStarFinder(
    # these defaults extracted from the
    # romancal SourceDetectionStep
    fwhm=1.0,
    threshold=0.0,
    sharplo=0.0,
    sharphi=1.0,
    roundlo=-1.0,
    roundhi=1.0,
    peakmax=None,
)


def create_gridded_psf_model(
    path_prefix,
    filt,
    detector,
    oversample=11,
    fov_pixels=9,
    sqrt_n_psfs=2,
    overwrite=False,
    buffer_pixels=100,
    instrument_options=None,
    logging_level=None,
):
    """
    Compute a gridded PSF model for one SCA via
    `~webbpsf.gridded_library.CreatePSFLibrary`.

    Parameters
    ----------
    path_prefix : str or Path-like
        Prefix to the output file path for the gridded PSF model
        FITS file. A suffix denoting the detector name will
        be appended.
    filt : str
        Filter name, starting with "F". For example: `"F184"`.
    detector : str
        Computed gridded PSF model for this SCA.
        Examples include: `"SCA01"` or `"SCA18"`.
    oversample : int, optional
        Oversample factor, default is 11. See WebbPSF docs for details [1]_.
        Choosing an odd number makes the pixel convolution more accurate.
    fov_pixels : int, optional
        Field of view width [pixels]. Default is 12.
        See WebbPSF docs for details [1]_.
    sqrt_n_psfs : int, optional
        Square root of the number of PSFs to calculate, distributed uniformly
        across the detector. Default is 4.
    overwrite : bool, optional
        Overwrite output file if one already exists. Default is False.
    buffer_pixels : int, optional
        Calculate a grid of PSFs distributed uniformly across the detector
        at least ``buffer_pixels`` away from the detector edges. Default is 100.
    instrument_options : dict, optional
        Instrument configuration options passed to WebbPSF.
        For example, WebbPSF assumes Roman pointing jitter consistent with
        mission specs by default, but this can be turned off with:
        ``{'jitter': None, 'jitter_sigma': 0}``.
    logging_level : str, optional
        Set logging level by name if not `None`, otherwise inherit from
        the romancal logger.

    Returns
    -------
    gridmodel : `photutils.psf.GriddedPSFModel`
        Gridded PSF model evaluated at several locations on one SCA.
    model_psf_centroids : list of tuples
        Pixel locations of the PSF models calculated for ``gridmodel``.

    References
    ----------
    .. [1] `WebbPSF documentation for `webbpsf.JWInstrument.calc_psf`
       <https://webbpsf.readthedocs.io/en/latest/api/webbpsf.JWInstrument.html#webbpsf.JWInstrument.calc_psf>`_

    """
    if int(sqrt_n_psfs) != sqrt_n_psfs:
        raise ValueError(f"`sqrt_n_psfs` must be an integer, got {sqrt_n_psfs}.")
    n_psfs = int(sqrt_n_psfs) ** 2

    # webbpsf appends "_sca??.fits" to the requested path:
    expected_output_path = f"{path_prefix}_{detector.lower()}.fits"

    # Choose pixel boundaries for the grid of PSFs:
    start_pix = 0
    stop_pix = 4096

    # Choose locations on detector for each PSF:
    pixel_range = np.linspace(
        start_pix + buffer_pixels, stop_pix - buffer_pixels, int(sqrt_n_psfs)
    )

    # generate PSFs over a grid of detector positions [pix]
    model_psf_centroids = [(int(x), int(y)) for y in pixel_range for x in pixel_range]

    if not os.path.exists(expected_output_path) or overwrite:
        if logging_level is None:
            # pass along logging level from __name__'s logger to WebbPSF:
            logging_level = logging.getLevelName(log.level)

        # set the WebbPSF logging level (similar to webbpsf.utils.setup_logging):
        conf.logging_level = logging_level
        restart_logging(verbose=False)

        wfi = webbpsf.roman.WFI()
        wfi.filter = filt

        if instrument_options is not None:
            wfi.options.update(instrument_options)

        # Initialize the PSF library
        inst = gridded_library.CreatePSFLibrary(
            instrument=wfi,
            filter_name=filt,
            detectors=detector.upper(),
            num_psfs=n_psfs,
            oversample=oversample,
            fov_pixels=fov_pixels,
            add_distortion=False,
            crop_psf=False,
            save=True,
            filename=path_prefix,
            overwrite=overwrite,
            verbose=False,
        )

        inst.location_list = model_psf_centroids

        # Create the PSF grid:
        gridmodel = inst.create_grid()

    elif os.path.exists(expected_output_path):
        logging.log(
            logging.INFO,
            f"Loading existing gridded PSF model from {expected_output_path}",
        )
        psf_model = CCDData.read(expected_output_path, unit=u.electron / u.s, ext=0)
        # the FITS file saved by webbpsf gets loaded with pixel
        # axes flipped in both dimensions, so flip them back after loading:
        psf_model.data = psf_model.data[::-1, ::-1]
        psf_model.meta = dict(psf_model.meta)
        psf_model.meta["oversampling"] = oversample
        psf_model.meta["grid_xypos"] = np.array(
            [list(tup)[::-1] for tup in model_psf_centroids]
        )
        gridmodel = GriddedPSFModel(psf_model)

    return gridmodel, model_psf_centroids


def fit_psf_to_image_model(
    image_model=None,
    data=None,
    error=None,
    dq=None,
    photometry_cls=PSFPhotometry,
    psf_model=None,
    grouper=None,
    fitter=None,
    localbkg_estimator=None,
    finder=None,
    x_init=None,
    y_init=None,
    progress_bar=False,
    error_lower_limit=None,
    fit_shape=(15, 15),
    exclude_out_of_bounds=True,
):
    """
    Fit PSF models to an ``ImageModel``.

    Parameters
    ----------
    image_model : `roman_datamodels.datamodels.ImageModel`
        Image datamodel. If ``image_model`` is supplied,
        ``data,error,dq`` should be `None`.
    data : `astropy.units.Quantity`
        Fit a PSF model to the rate image ``data``.
        If ``data,error,dq`` are supplied, ``image_model`` should be `None`.
    error : `astropy.units.Quantity`
        Uncertainties on fluxes in ``data``. Should be `None` if
        ``image_model`` is supplied.
    dq : `numpy.ndarray`
        Data quality bitmask for ``data``. Should be `None` if
        ``image_model`` is supplied.
    photometry_cls : {`photutils.psf.PSFPhotometry`,
            `photutils.psf.IterativePSFPhotometry`}
        Choose a photutils PSF photometry technique (default or iterative).
    psf_model : `astropy.modeling.Fittable2DModel`
        The 2D PSF model to fit to the rate image. Usually this model is an instance
        of `photutils.psf.GriddedPSFModel`.
    grouper : `photutils.psf.SourceGrouper`
        Specifies rules for attempting joint fits of multiple PSFs when
         there are nearby sources at small separations.
    fitter : `astropy.modeling.fitting.Fitter`, optional
        Modeling class which optimizes the PSF fit.
        Default is `astropy.modeling.fitting.LevMarLSQFitter(calc_uncertainties=True)`.
    localbkg_estimator : `photutils.background.LocalBackground`, optional
        Specifies inner and outer radii for computing flux background near
        a source. Default has ``inner_radius=10, outer_radius=30``.
    finder : subclass of `photutils.detection.StarFinderBase`, optional
        When ``photutils_cls`` is `photutils.psf.IterativePSFPhotometry`, the
        ``finder`` is called to determine if sources remain in the rate image
        after one PSF model is fit to the observations and removed.
        Default was extracted from the `DAOStarFinder` call in the
        Source Detection step.
    x_init : `numpy.ndarray`, optional
        Initial guesses for the ``x`` pixel coordinates of each source to fit.
    y_init : `numpy.ndarray`, optional
        Initial guesses for the ``y`` pixel coordinates of each source to fit.
    progress_bar : bool, optional
        Render a progress bar via photutils. Default is False.
    error_lower_limit : `astropy.units.Quantity`, optional
        Since some synthetic images may have bright sources with very
        small statistical uncertainties, the ``error`` can be clipped at
        ``error_lower_limit`` to prevent over-confident fits.
    fit_shape : int, or tuple of length 2, optional
        Rectangular shape around the center of a star that will
        be used to define the PSF-fitting data. See docs for
        `photutils.psf.PSFPhotometry` for details. Default is ``(16, 16)``.
    exclude_out_of_bounds : bool, optional
        If `True`, do not attempt to fit stars which have initial centroids
        that fall outside the pixel limits of the SCA. Default is False.

    Returns
    -------
    results_table : `astropy.table.QTable`
        PSF photometry results.
    photometry : instance of class ``photutils_cls``
        PSF photometry instance with configuration settings and results.

    """
    if grouper is None:
        # minimum separation before sources are fit simultaneously:
        grouper = SourceGrouper(min_separation=5)  # [pix]

    if fitter is None:
        fitter = LevMarLSQFitter(calc_uncertainties=True)

    # the iterative PSF method requires a finder:
    psf_photometry_kwargs = {}
    if photometry_cls is IterativePSFPhotometry or (x_init is None and y_init is None):
        if finder is None:
            finder = default_finder
        psf_photometry_kwargs["finder"] = finder

    if localbkg_estimator is not None:
        localbkg_estimator = LocalBackground(
            inner_radius=10,  # [pix]
            outer_radius=30,  # [pix]
        )

    photometry = photometry_cls(
        grouper=grouper,
        localbkg_estimator=localbkg_estimator,
        psf_model=psf_model,
        fitter=fitter,
        fit_shape=fit_shape,
        aperture_radius=fit_shape[0],
        progress_bar=progress_bar,
        **psf_photometry_kwargs,
    )

    if x_init is not None and y_init is not None:
        guesses = Table(np.column_stack([x_init, y_init]), names=["x_init", "y_init"])
    else:
        guesses = None

    if image_model is None:
        if data is None and error is None:
            raise ValueError(
                "PSF fitting requires either an ImageModel, "
                "or arrays for the data and error."
            )

    if dq is None:
        if image_model is not None:
            mask = dq_to_boolean_mask(image_model)
        else:
            mask = None
    else:
        mask = dq_to_boolean_mask(dq)

    if data is None and image_model is not None:
        data = image_model.data

    if error is None and image_model is not None:
        error = image_model.err

    if error_lower_limit is not None:
        # option to enforce a lower limit on the flux uncertainties
        error = np.clip(error, error_lower_limit, None)

    # we also mask non-finite values in the data and error arrays:
    non_finite = ~np.isfinite(data) | ~np.isfinite(error)

    if exclude_out_of_bounds and guesses is not None:
        # don't attempt to fit PSFs for objects with initial centroids
        # outside the detector boundaries:
        init_centroid_in_range = (
            (guesses["x_init"] > 0)
            & (guesses["x_init"] < data.shape[0])
            & (guesses["y_init"] > 0)
            & (guesses["y_init"] < data.shape[1])
        )
        guesses = guesses[init_centroid_in_range]

    # fit the model PSF to the data:
    results_table = photometry(
        data=data, error=error, init_params=guesses, mask=mask | non_finite
    )

    # results are stored on the PSFPhotometry instance:
    return results_table, photometry


def dq_to_boolean_mask(image_model_or_dq, ignore_flags=0, flag_map_name="ROMAN_DQ"):
    """
    Convert a DQ bitmask to a boolean mask. Useful for photutils methods.

    Parameters
    ----------
    image_model_or_dq : `roman_datamodels.datamodels.ImageModel` or `numpy.ndarray`
        ImageModel containing the DQ bitmask to convert to a boolean mask,
        or the DQ bitmask itself.
    ignore_flags : int, str, list, None (default = 0)
        See docs for `astropy.nddata.bitmask.extend_bit_flag_map`
    flag_map_name : str
        Name for the bitmask flag map in the astropy bitmask registry

    Returns
    -------
    mask : `numpy.ndarray`
        Boolean mask
    """

    if isinstance(image_model_or_dq, ImageModel):
        dq = image_model_or_dq.dq
    else:
        dq = image_model_or_dq

    # add the Roman DQ flags to the astropy bitmask registry:
    dq_flag_map = dict(roman_dq_flag_map)
    dq_flag_map.pop("GOOD")

    bitmask.extend_bit_flag_map(flag_map_name, **dq_flag_map)

    # convert the bitmask to a boolean mask:
    mask = bitmask.bitfield_to_boolean_mask(dq, ignore_flags=ignore_flags)
    return mask.astype(bool)
