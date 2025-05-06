"""
Module to calculate PSF photometry.
"""

import logging

import astropy.units as u
import numpy as np
import stpsf
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.table import Table
from astropy.utils import lazyproperty
from photutils.background import LocalBackground
from photutils.detection import DAOStarFinder
from photutils.psf import IterativePSFPhotometry, PSFPhotometry, SourceGrouper
from stpsf import gridded_library

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_gridded_psf_model(
    filt,
    detector,
    oversample=11,
    fov_pixels=9,
    sqrt_n_psfs=2,
    buffer_pixels=100,
    instrument_options=None,
):
    """
    Compute a gridded PSF model for one SCA via
    `~stpsf.gridded_library.CreatePSFLibrary`.

    Parameters
    ----------
    filt : str
        Filter name, starting with "F". For example: `"F184"`.
    detector : str
        Computed gridded PSF model for this SCA.
        Examples include: `"SCA01"` or `"SCA18"`.
    oversample : int, optional
        Oversample factor, default is 11. See STPSF docs for details [1]_.
        Choosing an odd number makes the pixel convolution more accurate.
    fov_pixels : int, optional
        Field of view width [pixels]. Default is 12.
        See STPSF docs for details [1]_.
    sqrt_n_psfs : int, optional
        Square root of the number of PSFs to calculate, distributed uniformly
        across the detector. Default is 4.
    buffer_pixels : int, optional
        Calculate a grid of PSFs distributed uniformly across the detector
        at least ``buffer_pixels`` away from the detector edges. Default is 100.
    instrument_options : dict, optional
        Instrument configuration options passed to STPSF.
        For example, STPSF assumes Roman pointing jitter consistent with
        mission specs by default, but this can be turned off with:
        ``{'jitter': None, 'jitter_sigma': 0}``.

    Returns
    -------
    gridmodel : `photutils.psf.GriddedPSFModel`
        Gridded PSF model evaluated at several locations on one SCA.
    model_psf_centroids : list of tuples
        Pixel locations of the PSF models calculated for ``gridmodel``.

    References
    ----------
    .. [1] `STPSF documentation for `stpsf.JWInstrument.calc_psf`
       <https://stpsf.readthedocs.io/en/latest/api/stpsf.JWInstrument.html#stpsf.JWInstrument.calc_psf>`_

    """
    if int(sqrt_n_psfs) != sqrt_n_psfs:
        raise ValueError(f"`sqrt_n_psfs` must be an integer, got {sqrt_n_psfs}.")
    n_psfs = int(sqrt_n_psfs) ** 2

    # Choose pixel boundaries for the grid of PSFs:
    start_pix = 0
    stop_pix = 4096

    # Choose locations on detector for each PSF:
    if sqrt_n_psfs != 1:
        pixel_range = np.linspace(
            start_pix + buffer_pixels, stop_pix - buffer_pixels, int(sqrt_n_psfs)
        )
    else:
        pixel_range = [(start_pix + stop_pix) / 2]

    # generate PSFs over a grid of detector positions [pix]
    model_psf_centroids = [(int(x), int(y)) for y in pixel_range for x in pixel_range]

    wfi = stpsf.roman.WFI()
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
        save=False,
        verbose=False,
    )

    inst.location_list = model_psf_centroids

    # Create the PSF grid:
    gridmodel = inst.create_grid()

    return gridmodel, model_psf_centroids


def fit_psf_to_image_model(
    image_model=None,
    data=None,
    error=None,
    mask=None,
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
        ``data,error`` should be `None`.
    data : `astropy.units.Quantity`
        Fit a PSF model to the rate image ``data``.
        If ``data,error`` are supplied, ``image_model`` should be `None`.
    error : `astropy.units.Quantity`
        Uncertainties on fluxes in ``data``. Should be `None` if
        ``image_model`` is supplied.
    mask : 2D bool `numpy.ndarray`, optional
        Mask to apply to the data. Default is `None`.
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
            # these defaults extracted from the
            # romancal SourceDetectionStep
            finder = DAOStarFinder(
                fwhm=1.0,
                threshold=0.0,
                sharplo=0.0,
                sharphi=1.0,
                roundlo=-1.0,
                roundhi=1.0,
                peakmax=None,
            )

        psf_photometry_kwargs["finder"] = finder

    if localbkg_estimator is None:
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

    if data is None and image_model is not None:
        data = image_model.data

    if error is None and image_model is not None:
        error = image_model.err

    if error_lower_limit is not None:
        # option to enforce a lower limit on the flux uncertainties
        error = np.clip(error, error_lower_limit, None)

    if exclude_out_of_bounds and guesses is not None:
        # don't attempt to fit PSFs for objects with initial centroids
        # outside the detector boundaries:
        init_centroid_in_range = (
            (guesses["x_init"] > 0)
            & (guesses["x_init"] < data.shape[1])
            & (guesses["y_init"] > 0)
            & (guesses["y_init"] < data.shape[0])
        )
        guesses = guesses[init_centroid_in_range]

    # fit the model PSF to the data:
    results_table = photometry(data=data, error=error, init_params=guesses, mask=mask)

    # results are stored on the PSFPhotometry instance:
    return results_table, photometry


class PSFCatalog:
    """
    Class to calculate PSF photometry.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        The input data model. The image data is assumed to be background
        subtracted.

    xypos : `numpy.ndarray`
        Pixel coordinates of sources to fit. The shape of this array should
        be (N, 2), where N is the number of sources to fit.

    mask : 2D `~numpy.ndarray` or `None`, optional
        A 2D boolean mask image with the same shape as the input data.
        This mask is used for PSF photometry. The mask should be the
        same one used to create the segmentation image.
    """

    def __init__(self, model, xypos, mask=None):
        self.model = model
        self.xypos = xypos
        self.mask = mask

        self.names = list(self._name_map.values())
        self.names.extend(["ra_psf", "dec_psf", "ra_psf_err", "dec_psf_err"])

        self.calc_psf_photometry()

    @lazyproperty
    def psf_model(self):
        """
        A gridded PSF model based on instrument and detector
        information.

        The `~photutils.psf.GriddedPSF` model is created using the
        STPSF library.
        """
        log.info("Constructing a gridded PSF model")
        if hasattr(self.model.meta, "instrument"):
            # ImageModel (L2 datamodel)
            filt = self.model.meta.instrument.optical_element
            detector = self.model.meta.instrument.detector.replace("WFI", "SCA")
        else:
            # MosaicModel (L3 datamodel)
            filt = self.model.meta.basic.optical_element
            detector = "SCA02"

        gridded_psf_model, _ = create_gridded_psf_model(
            filt=filt,
            detector=detector,
        )

        return gridded_psf_model

    @lazyproperty
    def _name_map(self):
        """
        Mapping of photutils column names to the output catalog names.
        """
        name_map = {}
        name_map["x_fit"] = "x_psf"
        name_map["x_err"] = "x_psf_err"
        name_map["y_fit"] = "y_psf"
        name_map["y_err"] = "y_psf_err"
        name_map["flux_fit"] = "psf_flux"
        name_map["flux_err"] = "psf_flux_err"
        name_map["qfit"] = "psf_gof"
        name_map["flags"] = "psf_flags"
        return name_map

    def calc_psf_photometry(self):
        """
        Perform PSF photometry by fitting PSF models to detected
        sources.
        """
        xinit, yinit = np.transpose(self.xypos)
        psf_photometry_table, _ = fit_psf_to_image_model(
            image_model=self.model,
            mask=self.mask,
            psf_model=self.psf_model,
            x_init=xinit,
            y_init=yinit,
            exclude_out_of_bounds=True,
        )

        # set these columns as attributes of this instance
        for old_name, new_name in self._name_map.items():
            value = psf_photometry_table[old_name]

            # change the photutils dtypes
            if np.issubdtype(value.dtype, np.integer):
                value = value.astype(np.int32)
            elif np.issubdtype(value.dtype, np.floating):
                value = value.astype(np.float32)

            # handle any unit conversions
            if new_name in ("x_psf", "y_psf", "x_psf_err", "y_psf_err"):
                value *= u.pix

            setattr(self, new_name, value)

    @lazyproperty
    def _sky_psf(self):
        """
        The pixel coordinates of the fitted-PSF pixel position.
        """
        return self.model.meta.wcs.pixel_to_world(self.x_psf, self.y_psf)

    @lazyproperty
    def ra_psf(self):
        """
        Right Ascension of the fitted-PSF pixel position.
        """
        return self._sky_psf.ra

    @lazyproperty
    def dec_psf(self):
        """
        Declination of the fitted-PSF pixel position.
        """
        return self._sky_psf.dec

    @lazyproperty
    def ra_psf_err(self):
        """
        The uncertainty in ra_psf.
        """
        return np.zeros(self.ra_psf.shape, dtype=np.float32) * u.arcsec

    @lazyproperty
    def dec_psf_err(self):
        """
        The uncertainty in dec_psf.
        """
        return np.zeros(self.dec_psf.shape, dtype=np.float32) * u.arcsec
