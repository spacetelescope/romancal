"""
Module to calculate PSF photometry.
"""

import logging
from collections import OrderedDict

import astropy.units as u
import numpy as np
from astropy.convolution import Box2DKernel, convolve
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
from astropy.table import Table
from astropy.utils import lazyproperty
from numpy import fft
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


def _cart_to_polar(image, oversample=2, order=3):
    """Convert an image from cartesian to polar coordinates.

    This makes a dense grid of points in polar coordinates and
    interpolates the values of image at those locations.

    Parameters
    ----------
    image : np.ndarray
        The image to convert to polar coordinates

    oversample : int
        The oversampling factor

    order : int
        The order of the interpolation used to extract
        points from the image

    Returns
    -------
    imagepolar, rr, tt
    imagepolar : np.ndarray
        the image in polar coordinates
    rr : np.ndarray
        the radii in pixel coordinates corresponding to
        the first dimension of imagepolar
    tt : np.ndarray
        the angles in radians corresponding to the
        second dimension of imagepolar
    """
    ntheta = 4 * image.shape[0] * oversample
    nrad = image.shape[0] * oversample
    szo2 = image.shape[0] // 2
    rr = np.linspace(0, szo2 * np.sqrt(2), nrad)
    tt = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    xx = rr[:, None] * np.cos(tt[None, :]) + szo2
    yy = rr[:, None] * np.sin(tt[None, :]) + szo2
    imagepolar = map_coordinates(
        image,
        [yy, xx],
        order=order,
        mode="constant",
        output=image.dtype,
        cval=np.nan,
    )
    return imagepolar, rr, tt


def _azimuthally_average_via_fft(data, pixel_scale_ratio=1.0):
    """Average a function azimuthally via FFT.

    Assumes that the data to be averaged is roughly circularly symmetric about
    the center of the image.  If not, large phase changes from pixel to pixel
    in the FFT will ruin the averaging in the FFT.

    Smoothing via the FFT has advantages over direct smoothing in real space
    when the image to be smoothed has nice properties in FFT space.
    For example, for diffraction-limited PSF stamps, the FFT goes to
    zero beyond some radius.  This averaging will preserve that behavior,
    which can otherwise be lost.

    Parameters
    ----------
    data : np.ndarray
        the image to be averaged

    pixel_scale_ratio : int
        the desired output pixel scale ratio divided by the input
        pixel scale ratio

    Returns
    -------
    image : np.ndarray
        the image after azimuthal averaging
    """
    if data.shape[0] != data.shape[1]:
        raise ValueError("image must be square")

    imfft = fft.fftshift(fft.fft2(fft.ifftshift(data)))
    polar, rr, tt = _cart_to_polar(imfft)
    kk = rr / data.shape[0]

    azimuthal_average = np.real(np.nanmean(polar, axis=1))
    npts = 2 * (int(np.ceil(data.shape[0] / pixel_scale_ratio)) // 2) + 1

    newk = fft.fftshift(fft.fftfreq(npts, d=pixel_scale_ratio))
    newkgrid = np.meshgrid(newk, newk)
    newkgrid_r = np.hypot(newkgrid[0], newkgrid[1])
    f_interp = np.interp(newkgrid_r.ravel(), kk, azimuthal_average, left=0, right=0)
    f_interp = f_interp.reshape(npts, npts)
    psf_new = np.real(fft.fftshift(fft.ifft2(fft.ifftshift(f_interp))))
    return psf_new


def _downsample_by_interpolation(image, downsample):
    """Downsample an image by interpolating it, preserving the centering.

    This is conceptually similar to taking every nth pixel of the image,
    but is careful about keeping the image centered and.  This is important
    for PSFs, for example, where we want to keep the PSF precisely centered,
    including for cases where the shape of the image and the amount of
    downsampling don't align neatly.

    We do this via linear interpolation, finding the locations in the
    original image that we want in the final downsampled image, and
    linearly interpolating to get the output values at those locations.

    Parameters
    ----------
    image : np.ndarray
        The image to be downsampled
    downsample : int
        the amount to downsample

    Returns
    -------
    image_downsampled : np.ndarray
        The downsampled image
    """
    ny, nx = image.shape
    ny_low = int(np.ceil(ny / downsample))
    nx_low = int(np.ceil(nx / downsample))

    # physical center of high-res image
    cy = (ny - 1) / 2
    cx = (nx - 1) / 2

    # offsets to low-res pixel centers
    y_offsets = (np.arange(ny_low) - (ny_low - 1) / 2) * downsample
    x_offsets = (np.arange(nx_low) - (nx_low - 1) / 2) * downsample

    # actual coordinates in high-res image
    y = cy + y_offsets
    x = cx + x_offsets
    yy, xx = np.meshgrid(y, x, indexing="ij")

    # interpolate
    low_res = map_coordinates(image, [yy, xx], order=1, mode="nearest")

    low_res *= image.sum() / low_res.sum()
    return low_res


def create_convolution_kernel(
    input_psf, target_psf, min_fft_power_ratio=1e-5, downsample=None, size=None
):
    """Find convolution kernel which convolves input_psf to match target_psf.

    The nominal photutils matching kernel code does a straight ratio
    of the target and input PSFs in Fourier space.  This is a little
    fraught for our PSFs where they go to ~0 at high frequencies, and
    unless the window function is also exactly zero there, you can get
    huge amounts of power.

    We mitigate this by taking an approach analogous to Boucaud+2016,
    adding a regularizing term to the denominator of the FFT.  The
    scale of this term is set by min_fft_power_ratio; its amplitude
    will be the maximum of the power in the input PSF times this
    ratio.  Large values correspond to stronger regularization.

    Boucaud+2016 penalizes variation in the matching kernel, while the
    approach taken here simply penalizes large values in the matching
    kernel.  Other approaches like penalizing second or third
    derivatives of the matching kernel could be taken, but for the
    light regularization needed here the amplitude of the matching
    kernel was adequate.

    Parameters
    ----------
    input_psf : np.ndarray
        The input PSF which needs to be convolved

    target_psf : np.ndarray
        The target PSF which input_psf should match following convolution

    min_fft_power_ratio : float
        controls the scale of the regularization of the matching kernel in
        terms of the peak power of the input PSF's FFT.

    downsample : int
        amount to downsample the final kernels by; useful if the models
        are oversampled but an image with the original sampling is desired

    size : int
        the desired size of the final stamp

    Returns
    -------
    kernel : np.ndarray
        Convolution kernel taking input_psf to target_psf

    """
    input_psf = fft.ifftshift(input_psf.copy())
    target_psf = fft.ifftshift(target_psf.copy())

    input_psf /= input_psf.sum()
    target_psf /= target_psf.sum()

    input_fft = fft.fft2(input_psf)
    target_fft = fft.fft2(target_psf)

    input_power = np.abs(input_fft) ** 2
    max_power = np.max(input_power)
    conv_kernel_fft = (
        target_fft
        * np.conj(input_fft)
        / (input_power + min_fft_power_ratio * max_power)
    )

    kernel = np.real(fft.fftshift(fft.ifft2(conv_kernel_fft)))
    kernel = kernel / kernel.sum()
    if downsample is not None and downsample > 1:
        kernel = convolve(kernel, Box2DKernel(width=downsample), boundary="extend")
        kernel = _downsample_by_interpolation(kernel, downsample)
    if size is not None:
        return central_stamp(kernel, size).copy()
    return kernel


def get_gridded_psf_model(psf_ref_model, focus=0, spectral_type=1):
    """Function to generate gridded PSF model from psf reference file

    Compute a gridded PSF model for one SCA using the
    reference files in CRDS.
    The input reference files have 3 focus positions and this is using
    the in-focus images. There are also three spectral types that are
    available and this code uses the M5V spectal type.

    Parameters
    ----------
    psf_ref_model : roman_datamodels.datamodels.EpsfRefModel
        PSF reference data model
    focus : integer
        index of slice into focus dimension of reference model.  0 for in focus
    spectral_type : integer
        index of slice into spectral type dimension of reference model.  1 for
        G2V.

    Returns
    -------
    photutils.psf.GriddedPSFModel of desired slice of ePSF reference file
    """
    # Open the reference file data model
    # select the infocus images (0) and we have a selection of spectral types
    # A0V, G2V, and M6V, pick G2V (1)
    psf_images = psf_ref_model.psf[focus, spectral_type, :, :, :].copy()
    # get the central position of the cutouts in a list
    psf_positions_x = psf_ref_model.meta.pixel_x
    psf_positions_y = psf_ref_model.meta.pixel_y
    meta = OrderedDict()
    position_list = []
    for index in range(len(psf_positions_x)):
        position_list.append([psf_positions_x[index], psf_positions_y[index]])

    # integrate over the native pixel scale
    oversample = psf_ref_model.meta.oversample
    pixel_response_kernel = Box2DKernel(width=oversample)
    for i in range(psf_images.shape[0]):
        psf = psf_images[i, :, :]
        im = convolve(psf, pixel_response_kernel) * oversample**2
        psf_images[i, :, :] = im

    meta["grid_xypos"] = position_list
    meta["oversampling"] = oversample
    nd = NDData(psf_images, meta=meta)
    model = GriddedPSFModel(nd)

    return model


def render_stamp(x, y, grid, size, oversample=1):
    """Render a PSF model at a specified location.

    Equivalent to grid.evaluate(xloc, yloc, 1, x, y) for xloc and yloc
    a square grid around (x, y) correspending to a [size, size] region.
    This function is meant to provide a more user-friendly PSF stamp
    rendering tool.

    Parameters
    ----------
    x : float
        x location at which to render PSF
    y : float
        y location at which to render PSF
    grid : photutils.psf.GriddedPSFModel
        gridded model to use for evaluation

    Returns
    -------
    stamp : np.ndarray[size, size]
        image of PSF at (x, y)
    """
    yy, xx = np.mgrid[:size, :size] / oversample

    # for a 3x3 stamp centered at [1, 1], we want [0, 1, 2], so
    # xcen and ycen should be 0.
    # for a 4x4 stamp, it's not clear where the stamp should be centered,
    # so either convention is fine.
    xcen, ycen = x - (size // 2) / oversample, y - (size // 2) / oversample
    stamp = grid.evaluate(xx + xcen, yy + ycen, 1, x, y)
    return stamp


def central_stamp(im, size):
    """Extract the central region of an image.

    The image must be square and we extract a square stamp.  The parity
    of the size of im and size must be the same; if they are not, we add
    one to size so that its parity matches that of the size of im.  The
    motivation for this behavior is that it's not clear, for example,
    what it means to extract the 'central' 1x1 pixel region from a 2x2
    stamp.

    Parameters
    ----------
    im : np.ndarray
        the image
    size : int
        the number of pixels to extract

    Returns
    -------
    stamp : np.ndarray[size, size]
        the central pixels of im, possibly adjusted by 1 to account
        for parity differences
    """

    if im.shape[0] != im.shape[1]:
        raise ValueError("im must be square")
    if (im.shape[0] % 2) != (size % 2):
        size = size + 1
    parity = int((im.shape[0] % 2) == 0)
    center = im.shape[0] // 2
    sizeo2 = size // 2
    return im[
        center - sizeo2 : center + sizeo2 + 1 - parity,
        center - sizeo2 : center + sizeo2 + 1 - parity,
    ]


def _get_jitter_params(meta):
    """Get jitter parameters out of metadata.

    Parameters
    ----------
    meta : metadata dictionary to extract keywords from.

    Returns
    -------
    dictionary containing jitter_major, jitter_minor, and jitter_position_angle
    keywords, substituting with sane defaults when not available.
    """
    return dict(
        jitter_major=meta.get("jitter_major", 8),
        jitter_minor=meta.get("jitter_minor", 8),
        jitter_position_angle=meta.get("jitter_position_angle", 0),
    )


def _evaluate_gaussian_fft(param, shape, pixel_scale):
    """Construct the FFT of a real-space Gaussian in the Fourier domain.

    We want to evaluate a jitter kernel in the Fourier domain.  The
    parameters of that jitter kernel are described in the 'param' dict.
    The jitter kernel needs to be convolved with PSF stamps from a
    reference file; the shape of those stamps and their PSF scale need
    to be specified in the shape and pixel_scale parameters.

    The jitter parameter dictionary must contain jitter_major, jitter_minor,
    and jitter_position_angle attributes.  The units of the major and
    minor axes must be in mas, while the position angle must be in degrees.
    The position angle is the angle of the major axis relative to the +Y
    axis, with PA increasing as the major axis moves toward the +X axis.
    Major and minor axis sizes are expressed in units of RMS; i.e., the
    standard deviation of the corresponding Gaussian.

    Because Gaussians are real-valued, only the positive frequency
    components are returned.

    Parameters
    ----------
    param : dict
        jitter parameter dictionary from _get_jitter_params
    shape : tuple
        shape of PSF stamp
    pixel_scale : float
        pixel scale of PSF, arcsec

    Returns
    -------
    fft : np.ndarray
        analytically computed FFT of real-space Gaussian
    """
    # pixel scale must include oversampling factor
    # i.e., be something like the native scale divided by 4

    # need to shift theta -> np.pi / 2 - theta in order to go from
    # 'math' convention of +x axis, increasing toward +y axis
    # to 'astronomy' convention of +y axis, increasing toward
    # +x axis.
    theta = np.pi / 2 - np.radians(param["jitter_position_angle"])
    major = param["jitter_major"] / 1000
    minor = param["jitter_minor"] / 1000
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    # rfftfreq since this is real
    uu = np.fft.rfftfreq(shape[1], d=pixel_scale)
    vv = np.fft.fftfreq(shape[0], d=pixel_scale)
    ugrid, vgrid = np.meshgrid(uu, vv)
    covar = rot.dot(np.diag([major**2, minor**2])).dot(rot.T)
    uvgrid = np.array([ugrid.reshape(-1), vgrid.reshape(-1)])
    rr2 = np.einsum("ji,jk,ki->i", uvgrid, covar, uvgrid)
    fft = np.exp(-2 * np.pi**2 * rr2)
    fft = fft.reshape(ugrid.shape)
    return fft


def add_jitter(psf_ref_model, image_model, pixel_scale=0.11):
    """Add jitter to PSF stamps in a PSF reference file.

    The PSF reference model contains some amount of jitter.  The image
    model contains some potentially different amount of jitter.  This
    function needs to deconvolve the jitter in the reference model
    and then reconvolve with the jitter in the image, to obtain a model
    for the jitter in this image.

    This procedure is made more challenging by the fact that the jitter is
    small.  Typical values of the jitter are 8 mas compared to a typical
    PSF FWHM of ~100 mas.  The pixel scale in the reference files is 4x
    oversampled relative to the native scale of 0.11", and so is ~28 mas---
    much larger than the 8 mas jitter scale, so naive direct convolution
    will not work.

    We use the approach of Lang (2020) to do the convolutions in the Fourier
    domain to avoid this problem.  The key element there is to analytically
    evaluate the jitter kernel in Fourier space rather than to evaluate
    it in real space and then transform it.  If the PSF is well sampled,
    all fourier components corresponding to higher frequencies than the Nyquist
    frequency will be exactly zero.  Accordingly, if we evaluate the jitter
    kernel in the FFT domain analytically only over the frequencies where the
    PSF FFT is non-zero, we can multiply the two FFTs and transform back
    to real space without needing to think about the fact that the jitter
    kernel is undersampled in real space.

    Parameters
    ----------
    psf_ref_model : datamodels.EpsfRefModel
        ePSF reference models from CRDS
    image_model : datamodels.ImageModel
        image for which a PSF model is desired
    pixel_scale : float
        pixel scale of the image.  Because the jitter kernel is very small,
        this does not need to be terribly accurate; it's used to transform
        the real-space jitter kernel specification into the pixel-space
        kernels that are convolved with the PSF stamps.

    Returns
    -------
    convolved_stamps : np.ndarray
        PSF reference file stamps after convolution with jitter kernel
    """

    ref_jitter_params = _get_jitter_params(psf_ref_model.meta)
    image_jitter_params = _get_jitter_params(
        getattr(image_model.meta, "guide_star", dict())
    )
    oversample = psf_ref_model.meta.oversample
    shape = psf_ref_model.psf.shape[-2:]
    jit_ref_fft = _evaluate_gaussian_fft(
        ref_jitter_params, shape, pixel_scale / oversample
    )
    jit_img_fft = _evaluate_gaussian_fft(
        image_jitter_params, shape, pixel_scale / oversample
    )
    jitter_fft = jit_img_fft / jit_ref_fft
    stamp_fft = np.fft.rfft2(psf_ref_model.psf, axes=(-2, -1))
    convolved = np.fft.irfft2(stamp_fft * jitter_fft[None, ...], axes=(-2, -1), s=shape)
    log.info(
        f"Performed jitter convolution.  Image kernel: {tuple(image_jitter_params.values())} Reference kernel: {tuple(ref_jitter_params.values())}"
    )
    return convolved


def create_l3_psf_model(
    psf_ref_model,
    stamp_radius=10,
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
    detector_pixel_scale = 0.11  # Roman native scale
    # psf is now a stamp going 10 pixels out from the center of the PSF
    # at the native PSF scale, oversampled by a factor of oversample.

    # Smooth to account for the pixfrac used to create the L3 image.
    pixfrac_kernel = Box2DKernel(width=pixfrac * oversample)
    psf = convolve(psf, kernel=pixfrac_kernel)
    # Smooth to the image scale
    outscale_kernel = Box2DKernel(width=oversample * pixel_scale / detector_pixel_scale)
    psf = convolve(psf, kernel=outscale_kernel)
    # Azimuthally smooth the psf
    psf = _azimuthally_average_via_fft(
        psf, pixel_scale_ratio=pixel_scale / detector_pixel_scale
    )

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
