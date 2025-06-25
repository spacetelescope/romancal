"""
Module to calculate PSF photometry.
"""

import logging
import math
from collections import OrderedDict

import astropy.units as u
import numpy as np
import scipy.ndimage as ndimage
from astropy.convolution import Box2DKernel, convolve
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
from astropy.table import Table
from astropy.utils import lazyproperty
from photutils.background import LocalBackground
from photutils.detection import DAOStarFinder
from photutils.psf import (
    GriddedPSFModel,
    ImagePSF,
    IterativePSFPhotometry,
    PSFPhotometry,
    SourceGrouper,
)
from scipy.ndimage import map_coordinates

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def azimuthally_smooth(data, oversample=2, scaling=1.0, order=4):
    """Azimuthally smooth model

    The image is converted to polar coordinates via a 4th order spline interpolation.
    The image average is determined at each radius, and a final image is constructed by reprojecting this averaged
    image back into cartesian coordinates.

    Parameters
    ----------
    data : nd.array(size=(*, *))
        Data to be smoothed.

    oversample : int
        Oversampling of the data to improve fidelity of the conversions
        between cartesian and polar layout

    scaling : float
        Scale factor to apply to the result.

    order : int
        Order of the spline interpolation used for the coordinate transformations

    Returns
    -------
    smoothed : nd.array(size=(*, *))
        The azimuthally smoothed image
    """

    # Define cartesian->polar->cartesian conversion.
    # Significant difference is the polar back to cartesian uses
    # the original data's transformation parameters, to reproduce
    # an image of the correct dimensions.
    def cart_to_polar(model, oversample=2, order=4):
        ntheta = 4 * model.shape[0] * oversample
        nrad = model.shape[0] * oversample
        szo2 = model.shape[0] // 2
        rr = np.linspace(-5, szo2 * np.sqrt(2), nrad)
        tt = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
        xx = rr[:, None] * np.cos(tt[None, :]) + szo2
        yy = rr[:, None] * np.sin(tt[None, :]) + szo2
        modelpolar = map_coordinates(
            model, [xx, yy], order=order, mode="nearest", output=np.dtype("f4")
        )
        return modelpolar, rr, tt

    def polar_to_cart(model, rr, tt, scaling=1.0, order=4):
        szo = model.shape[0]
        szo = math.ceil(szo / scaling)
        szo2 = szo // 2
        xo, yo = np.mgrid[-szo2 : szo2 + 1, -szo2 : szo2 + 1] * scaling
        ro = np.sqrt(xo**2 + yo**2)
        to = np.arctan2(yo, xo) % (2 * np.pi)
        ro = np.interp(ro, rr, np.arange(len(rr)).astype("f4"))
        to = np.interp(to, tt, np.arange(len(tt)).astype("f4"))
        res = map_coordinates(
            model, [ro, to], order=order, mode="wrap", output=np.dtype("f4")
        )
        return res

    polar, rr, tt = cart_to_polar(data, oversample=oversample)
    polar_mean = np.mean(polar, axis=1, keepdims=True)
    smoothed = polar_to_cart(polar_mean, rr, tt, scaling=scaling)

    return smoothed


def get_psf_library(self):
    """Function to retrieve psf library from CRDS

    Compute a gridded PSF model for one SCA using the
    reference files in CRDS.
    The input reference files have 3 focus positions and this is using
    the in-focus images. There are also three spectral types that are
    available and this code uses the M5V spectal type.
    """
    # Open the reference file data model
    # select the infocus images (0) and we have a selection of spectral types
    # A0V, G2V, and M6V, pick M5V (2)
    focus = 0
    spectral_type = 2
    jitter_value = 1.1
    psf_images = self.psf_ref_model.psf[focus, spectral_type, :, :, :]
    psf_images = ndimage.gaussian_filter(psf_images, sigma=jitter_value)
    # get the central position of the cutouts in a list
    psf_positions_x = self.psf_ref_model.meta.pixel_x.data.data
    psf_positions_y = self.psf_ref_model.meta.pixel_y.data.data
    meta = OrderedDict()
    position_list = []
    for index in range(len(psf_positions_x)):
        position_list.append([psf_positions_x[index], psf_positions_y[index]])

    meta["grid_xypos"] = position_list
    meta["oversampling"] = self.psf_ref_model.meta.oversample
    nd = NDData(psf_images, meta=meta)
    model = GriddedPSFModel(nd)

    return model


def create_l3_psf_model(
    filt,
    detector="SCA02",
    pixfrac=1.0,
    pixel_scale=0.11,
    oversample=11,
    fov_pixels=9,
    instrument_options=None,
):
    """
    Compute a PSF model via `~stpsf.calc_psf`.

    L3 data is an amalgamation of numerous exposures over numerous SCA's.
    This algorithm does not attempt to merge specific PSF profiles for each
    SCA that contributes to each pixel. Instead, a simplified version is implemented
    as described.

        - Base PSF for the given detector, is created. This base has the
          0.11 arcsec pixel scale convolved in.
        - PSF is further convolved with the drizzlepac `pixfrac` scale
        - PSF is further convolved with the images actual pixel scale.
        - The PSF is then azimuthally averaged and resampled at the L3 pixel scale.

    Parameters
    ----------
    filt : str
        Filter name, starting with "F". For example: `"F184"`.
    detector : str
        Computed gridded PSF model for this SCA.
        Examples include: `"SCA01"` or `"SCA18"`.
    pixfrac : float
        drizzlepac pixel fraction used.
    pixel_scale : float
        L3 image pixel scale in arcsec.
        Often similar to the default detector scale of 0.11 arcsec.
    oversample : int, optional
        Oversample factor, default is 11. See STPSF docs for details [1]_.
        Choosing an odd number makes the pixel convolution more accurate.
    fov_pixels : int, optional
        Field of view width [pixels]. Default is 12.
        See STPSF docs for details [1]_.
    instrument_options : dict, optional
        Instrument configuration options passed to STPSF.
        For example, STPSF assumes Roman pointing jitter consistent with
        mission specs by default, but this can be turned off with:
        ``{'jitter': None, 'jitter_sigma': 0}``.

    Returns
    -------
    psf_model : `photutils.psf.ImagePSF`
        PSF model.

    References
    ----------
    .. [1] `STPSF documentation for `stpsf.JWInstrument.calc_psf`
       <https://stpsf.readthedocs.io/en/latest/api/stpsf.JWInstrument.html#stpsf.JWInstrument.calc_psf>`_

    """

    # Create base PSF.
    wfi = WFI()
    wfi.detector = detector
    wfi.filter = filt
    wfi_psf = wfi.calc_psf(
        fov_pixels=fov_pixels,
        oversample=oversample,
        add_distortion=False,
        crop_psf=False,
    )
    psf = wfi_psf[0].data
    detector_pixel_scale = wfi_psf[1].header["PIXELSCL"]

    # Pixel response
    pixel_response_kernal = Box2DKernel(width=oversample)
    psf = convolve(psf, pixel_response_kernal)

    # Smooth to account for the pixfrac used to create the L3 image.
    pixfrac_kernel = Box2DKernel(width=pixfrac * oversample)
    psf = convolve(psf, kernel=pixfrac_kernel)

    # Smooth to the image scale
    outscale_kernel = Box2DKernel(width=oversample * pixel_scale / detector_pixel_scale)
    psf = convolve(psf, kernel=outscale_kernel)

    # Azimuthally smooth the psf
    psf = azimuthally_smooth(psf, scaling=pixel_scale / detector_pixel_scale)

    # Create the PSF model.
    x_0, y_0 = psf.shape
    x_0 = (x_0 - 1) / 2.0 / oversample
    y_0 = (y_0 - 1) / 2.0 / oversample
    psf_model = ImagePSF(psf, x_0=x_0, y_0=y_0, oversampling=oversample)

    return psf_model


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

    def __init__(self, model, psf_ref_model, xypos, mask=None):
        self.model = model
        self.psf_ref_model = psf_ref_model
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
        """
        if hasattr(self.model.meta, "instrument"):
            # ImageModel (L2 datamodel)
            gridded_psf_model = get_psf_library(self)
        else:
            raise NotImplementedError("Need to update l3 psf to use crds")
            # MosaicModel (L3 datamodel)
            filt = self.model.meta.basic.optical_element
            psf_model = create_l3_psf_model(
                filt=filt,
                pixel_scale=self.model.meta.wcsinfo.pixel_scale
                * 3600.0,  # wcsinfo is in degrees. Need arcsec
                pixfrac=self.model.meta.resample.pixfrac,
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
