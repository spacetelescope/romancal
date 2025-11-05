"""
Module to calculate PSF photometry.
"""

import logging
import math
from collections import OrderedDict

import astropy.units as u
import numpy as np
from numpy import fft
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
from photutils.psf import matching

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)



def cart_to_polar(image, oversample=2, order=3):
    """Convert an image from cartesian to polar coordinates.

    Parameters
    ----------
    image : np.ndarray
        The image to convert to polar coordinates

    oversample : int
        The oversampling factor 
    """
    ntheta = 4 * image.shape[0] * oversample
    nrad = image.shape[0] * oversample
    szo2 = image.shape[0] // 2
    rr = np.linspace(0, szo2 * np.sqrt(2), nrad)
    tt = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    xx = rr[:, None] * np.cos(tt[None, :]) + szo2
    yy = rr[:, None] * np.sin(tt[None, :]) + szo2
    imagepolar = map_coordinates(
        image, [yy, xx], order=order, mode="constant",
        output=image.dtype, cval=np.nan,
    )
    return imagepolar, rr, tt


def azimuthally_average_via_fft(data, oversample=2, pixel_scale_ratio=1.0):
    """Average a function azimuthally via FFT.

    Assumes that the data to be smooth is roughly circularly symmetric about
    the center of the image.  If not, large phase changes from pixel to pixel
    in the FFT will ruin the averaging in the FFT.

    Smoothing via the FFT rather than directly in real space is nice
    when the image to be smoothed has nice properties in FFT space.
    For example, for diffraction-limited PSF stamps, the FFT goes to
    zero beyond some radius.  This averaging will preserve that behavior,
    which can be lost in direct-space averaging approaches.
    """
    assert data.shape[0] == data.shape[1]

    imfft = fft.fftshift(fft.fft2(fft.ifftshift(data)))
    polar, rr, tt = cart_to_polar(imfft)
    kk = rr / data.shape[0]

    azimuthal_average = np.real(np.nanmean(polar, axis=1))
    npts = 2 * (int(np.ceil(data.shape[0] / pixel_scale_ratio)) // 2) + 1

    newk = fft.fftshift(fft.fftfreq(npts, d=pixel_scale_ratio))
    newkgrid = np.meshgrid(newk, newk)
    newkgrid_r = np.hypot(newkgrid[0], newkgrid[1])
    f_interp = np.interp(newkgrid_r.ravel(), kk, azimuthal_average, left=0, right=0)  # better interpolation needed?
    f_interp = f_interp.reshape(npts, npts)

    psf_new = np.real(fft.fftshift(fft.ifft2(fft.ifftshift(f_interp))))
    return psf_new


def create_convolution_kernel(input_psf, target_psf,
                                        min_fft_power_ratio=1e-3):
    """Find convolution kernel which convolves input_psf to match target_psf.

    The nominal photutils matching kernel code does a straight ratio of the target
    and input PSFs in Fourier space.  This is a little fraught for our PSFs
    where they go to ~0 at high frequencies, and unless the window function is
    also exactly zero there, you can get huge amounts of power.

    To mitigate this, we introduce a parameter that controls what FFT powers to
    set to zero due to being too small.  Values in the FFT smaller than this
    are zeroed.  If the input FFT is zero where the target FFT is non-zero,
    this function will emit an error message but will otherwise let those
    frequencies propagate through unchanged.  This will not produce a good matching
    PSF.

    Parameters
    ----------
    input_psf : np.ndarray
        The input PSF which needs to be convolved

    target_psf : np.ndarray
        The target PSF which input_psf should match following convolution

    min_fft_power_ratio : float
        zero power in the FFT when it is smaller than this ratio times its
        peak

    Returns
    -------
    kernel : np.ndarray
        Convolution kernel taking input_psf to target_psf
    """
    # if window is None:
    #     window = matching.TopHatWindow(0.25)

    input_psf = input_psf.copy()
    target_psf = target_psf.copy()
    
    input_psf /= input_psf.sum()
    target_psf /= target_psf.sum()

    input_fft = fft.fft2(fft.ifftshift(input_psf))
    target_fft = fft.fft2(fft.ifftshift(target_psf))

    target_power = np.abs(target_fft)
    input_power = np.abs(input_fft)
    mtarget = target_power > np.max(target_power) * min_fft_power_ratio
    minput = input_power > np.max(input_power) * min_fft_power_ratio
    target_fft[~mtarget] = 0
    input_fft[~minput] = 0

    if np.any((target_fft != 0) & (input_fft == 0)):
        log.error('Could not create good matching kernel!')
    
    
    ratio = target_fft / (input_fft + (input_fft == 0))

    kernel = fft.fftshift(np.real(fft.fftshift(fft.ifft2(ratio))))
    return kernel / kernel.sum()


    return matching.create_matching_kernel(
        input_psf.data, target_psf.data, window=window)


def get_gridded_psf_model(psf_ref_model):
    """Function to generate gridded PSF model from psf reference file

    Compute a gridded PSF model for one SCA using the
    reference files in CRDS.
    The input reference files have 3 focus positions and this is using
    the in-focus images. There are also three spectral types that are
    available and this code uses the M5V spectal type.
    """
    # Open the reference file data model
    # select the infocus images (0) and we have a selection of spectral types
    # A0V, G2V, and M6V, pick G2V (1)
    focus = 0
    spectral_type = 1
    psf_images = psf_ref_model.psf[focus, spectral_type, :, :, :].copy()
    # get the central position of the cutouts in a list
    psf_positions_x = psf_ref_model.meta.pixel_x.data.data
    psf_positions_y = psf_ref_model.meta.pixel_y.data.data
    meta = OrderedDict()
    position_list = []
    for index in range(len(psf_positions_x)):
        position_list.append([psf_positions_x[index], psf_positions_y[index]])

    # integrate over the native pixel scale
    oversample = psf_ref_model.meta.oversample
    pixel_response_kernel = Box2DKernel(width=oversample)
    for i in range(psf_images.shape[0]):
        psf = psf_images[i, :, :].copy()
        im = convolve(psf, pixel_response_kernel) * oversample**2
        psf_images[i, :, :] = im

    meta["grid_xypos"] = position_list
    meta["oversampling"] = oversample
    nd = NDData(psf_images, meta=meta)
    model = GriddedPSFModel(nd)

    return model


def create_l3_psf_model(
    psf_ref_model,
    stamp_radius = 10,
    pixfrac=1.0,
    pixel_scale=0.11,
    oversample=None,
):
    """
    Compute a PSF model for an L3 image.

    L3 data is an amalgamation of numerous exposures over numerous SCA's.
    This algorithm does not attempt to merge specific PSF profiles for each
    SCA that contributes to each pixel. Instead, a simplified version is implemented
    as described.

        - Base PSF for the given detector, is created via get_gridded_psf_model.
          This bas has the native 0.11 arcsec pixel scale already convolved in.
        - PSF is further convolved with the drizzlepac `pixfrac` scale
        - PSF is further convolved with the images actual pixel scale.
        - The PSF is then azimuthally averaged and resampled at the L3 pixel scale.

    Parameters
    ----------
    psf_ref_model : str
        PSF reference file model to use
    stamp_radius : int
        PSF model should extend this many native pixels from the center pixel
    pixfrac : float
        drizzlepac pixel fraction used.
    pixel_scale : float
        L3 image pixel scale in arcsec.
        Often similar to the default detector scale of 0.11 arcsec.
    oversample : int, optional
        Oversample factor, default uses gridpsf oversampling

    Returns
    -------
    psf_model : `photutils.psf.ImagePSF`
        PSF model.

    """
    gridpsf = get_gridded_psf_model(psf_ref_model)  # 361x361, 4x oversampled
    # this already includes integration over the native pixel scale
    center = 2044  # Roman SCA center pixel
    oversample = oversample if oversample is not None else gridpsf.oversampling[0]
    npix = stamp_radius * 2 * oversample + 1
    pts = np.linspace(-stamp_radius + center, stamp_radius + center, npix)
    xx, yy = np.meshgrid(pts, pts)
    psf = gridpsf.evaluate(xx, yy, 1, center, center)
    psf1 = psf.copy()
    detector_pixel_scale = 0.11  # Roman native scale
    # psf is now a stamp going 10 pixels out from the center of the PSF
    # at the native PSF scale, oversampled by a factor of oversample.

    # Smooth to account for the pixfrac used to create the L3 image.
    pixfrac_kernel = Box2DKernel(width=pixfrac * oversample)
    psf = convolve(psf, kernel=pixfrac_kernel)
    psf2 = psf.copy()
    # Smooth to the image scale
    outscale_kernel = Box2DKernel(width=oversample * pixel_scale / detector_pixel_scale)
    psf = convolve_box_tophat(psf, size=oversample * pixel_scale / detector_pixel_scale)
    psf3 = psf.copy()
    # Azimuthally smooth the psf
    psf = azimuthally_average_via_fft(psf, pixel_scale_ratio=pixel_scale / detector_pixel_scale)

    # Create the PSF model.
    x_0, y_0 = psf.shape
    x_0 = (x_0 - 1) / 2.0 / oversample
    y_0 = (y_0 - 1) / 2.0 / oversample
    psf_model = ImagePSF(psf, x_0=x_0, y_0=y_0, oversampling=oversample)
    import pdb
    pdb.set_trace()

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
            psf_model = get_gridded_psf_model(self.psf_ref_model)
        else:
            # MosaicModel (L3 datamodel)
            psf_model = create_l3_psf_model(
                self.psf_ref_model,
                pixel_scale=self.model.meta.wcsinfo.pixel_scale
                * 3600.0,  # wcsinfo is in degrees. Need arcsec
                pixfrac=self.model.meta.resample.pixfrac,
            )

        return psf_model

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
