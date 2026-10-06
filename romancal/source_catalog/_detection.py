"""
Module to detect sources using image segmentation.
"""

import logging
import math
import warnings

import numpy as np
import scipy.ndimage
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.utils.exceptions import AstropyUserWarning
from photutils.segmentation import SourceFinder, make_2dgaussian_kernel
from photutils.utils.exceptions import NoDetectionsWarning

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def convolve_data(data, kernel_fwhm, size=None, mask=None):
    """
    Convolve the background-subtracted model image with a Gaussian
    kernel.

    Parameters
    ----------
    data : 2D `numpy.ndarray`
        The input 2D array. The data is assumed to be
        background subtracted.
    kernel_fwhm : float
        The FWHM of the Gaussian kernel.
    size : int, optional
        The size of the kernel array. Default is 3 times the kernel
        FWHM.
    mask : 2D `numpy.ndarray`, optional
        A boolean mask with the same shape as the input data. True
        values indicate masked pixels.

    Returns
    -------
    convolved_data : `numpy.ndarray`
        The convolved data array.
    """
    size = math.ceil(kernel_fwhm * 3)
    size = size + 1 if size % 2 == 0 else size  # make size be odd
    kernel = make_2dgaussian_kernel(kernel_fwhm, size=size)  # normalized to 1

    with warnings.catch_warnings():
        # suppress warnings caused by large masked areas
        warnings.simplefilter("ignore", AstropyUserWarning)
        return convolve(data, kernel, mask=mask)


def make_segmentation_image(
    convolved_data,
    snr_threshold,
    n_pixels,
    bkg_rms,
    deblend=False,
    mask=None,
    deblend_contrast=1e-4,
):
    """
    Make a segmentation image from a model image.

    Parameters
    ----------
    convolved_data : 2D `numpy.ndarray`
        The background-subtracted convolved data array.
    snr_threshold : float
        The per-pixel signal-to-noise ratio threshold for detection.
    n_pixels : int
        The number of connected pixels required to define a source.
    bkg_rms : 2D `numpy.ndarray`
        The background RMS array.
    deblend : bool, optional
        Whether to deblend overlapping sources. Default is False.
    mask : 2D `numpy.ndarray`, optional
        A boolean mask with the same shape as the input data. True
        values indicate masked pixels.
    deblend_contrast : float, optional
        The fraction of the total source flux that a local peak must have to
        be deblended as a separate source.

    Returns
    -------
    segment_img : `SegmentationImage`
        The segmentation image.
    """
    with warnings.catch_warnings():
        # suppress NoDetectionsWarning from photutils
        warnings.filterwarnings("ignore", category=NoDetectionsWarning)

        finder = SourceFinder(n_pixels, deblend=deblend, contrast=deblend_contrast)
        threshold = snr_threshold * bkg_rms
        segment_img = finder(convolved_data, threshold, mask=mask)

        if segment_img is None:
            n_sources = 0
        else:
            n_sources = segment_img.n_labels
        log.info(f"Detected {n_sources} sources")

    return segment_img


def make_gaussian_kernel_1d(fwhm, size_factor=4, max_size=601):
    """
    Make a normalized 1D Gaussian kernel.

    The outer product of this kernel with itself is a normalized 2D
    circular Gaussian, so a 2D convolution with that Gaussian is two 1D
    convolutions with this one.  Each pixel is the mean of the Gaussian
    over a 10x oversampled grid, matching the photutils 2D kernels.

    Parameters
    ----------
    fwhm : float
        Full-width at half-maximum of the Gaussian, in pixels.
    size_factor : int, optional
        Kernel length as a multiple of ``fwhm``.
    max_size : int, optional
        Largest kernel length, in pixels.

    Returns
    -------
    kernel : 1D `numpy.ndarray`
        The kernel array, summing to 1.
    """
    size = min(math.ceil(size_factor * fwhm), max_size)
    size = size + 1 if size % 2 == 0 else size  # make size be odd
    kernel = Gaussian1DKernel(
        fwhm * gaussian_fwhm_to_sigma, x_size=size, mode="oversample", factor=10
    )
    kernel.normalize(mode="integral")
    return kernel.array


def make_gaussian_kernel(fwhm, size_factor=4, max_size=601):
    """
    Make a normalized 2D circular Gaussian kernel.

    Parameters
    ----------
    fwhm : float
        Full-width at half-maximum of the Gaussian, in pixels.
    size_factor : int, optional
        Kernel box size as a multiple of ``fwhm``.
    max_size : int, optional
        Largest kernel size, in pixels.

    Returns
    -------
    kernel : 2D `numpy.ndarray`
        The kernel array, summing to 1.
    """
    kernel = make_gaussian_kernel_1d(fwhm, size_factor=size_factor, max_size=max_size)
    return np.outer(kernel, kernel)


def separable_convolve(data, kernel, mask=None):
    """
    Convolve ``data`` with the outer product of ``kernel`` with itself.

    The convolution is done directly, as one 1D convolution along each
    axis, and is zero-padded at the boundary.  Unlike an FFT, a direct
    convolution is exactly zero far from any nonzero data, so masked
    regions and image edges are not filled with round-off from bright
    sources elsewhere in the image.

    This uses single precision to conserve memory and fills NaNs to zeros
    before convolution.

    Parameters
    ----------
    data : 2D `numpy.ndarray`
        The array to convolve.
    kernel : 1D `numpy.ndarray`
        The 1D convolution kernel, with odd length so that the result is
        centered.
    mask : 2D `numpy.ndarray`, optional
        Boolean mask; where True, data is zeroed before convolution.

    Returns
    -------
    convolved : 2D `numpy.ndarray`
        The convolved array, of the same shape as ``data``.
    """
    if mask is not None:
        data = np.where(mask, 0.0, data)
    data = np.nan_to_num(
        np.asarray(data, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0
    )
    kernel = np.asarray(kernel, dtype=np.float32)
    out = scipy.ndimage.convolve1d(data, kernel, axis=0, mode="constant")
    return scipy.ndimage.convolve1d(out, kernel, axis=1, mode="constant")


def ivw_convolve(data, wht, kernel, mask=None):
    """
    Inverse-variance-weighted convolution of data with uncertainties given a template.

    Construction of a significance image for a matched filter requires
    a signal and uncertainty image for each template.  This function computes
    both.  The ratio ``num / sqrt(denom2)`` returned by this function is the
    maximum-likelihood amplitude of a template kernel divided by its uncertainty,
    given some Gaussian noise.

    Masked pixels enter the convolution with zero weight, so the result is
    defined on masked pixels.

    Parameters
    ----------
    data : 2D `numpy.ndarray`
        Background-subtracted data.
    wht : 2D `numpy.ndarray`
        Inverse-variance weights; zero where masked.
    kernel : 1D `numpy.ndarray`
        The 1D factor of the separable template; see `separable_convolve`.
    mask : 2D `numpy.ndarray`, optional
        Boolean mask; True values force the data and weights to zero.

    Returns
    -------
    num : 2D `numpy.ndarray`
        ``conv(data * wht, kernel)``.
    denom2 : 2D `numpy.ndarray`
        ``conv(wht, kernel**2)``.
    """
    return (
        separable_convolve(data * wht, kernel, mask=mask),
        separable_convolve(wht, kernel**2, mask=mask),
    )


def snr_from_ivw(num, denom2):
    """
    Form a signal-to-noise ratio image from accumulated `ivw_convolve` arrays.

    Kept separate from `ivw_convolve` because the arrays must be summed
    over bands before the ratio is taken; templates, by contrast, are
    combined after it, by taking the maximum of their ratios.
    ``denom2`` is mathematically
    non-negative, but is clamped before the square root to guard against
    round-off.

    Parameters
    ----------
    num, denom2 : 2D `numpy.ndarray`
        The numerator and squared-denominator arrays, summed over bands.

    Returns
    -------
    snr : 2D `numpy.ndarray`
        ``num / sqrt(denom2)``, zero where there is no weight.
    """
    denom = np.sqrt(np.maximum(denom2, 0))
    return np.divide(num, denom, out=np.zeros_like(num), where=denom > 0)
