"""
Module to detect sources using image segmentation.
"""

import logging
import math
import warnings

import photutils
from astropy.convolution import convolve
from astropy.utils import minversion
from astropy.utils.exceptions import AstropyUserWarning
from photutils.segmentation import SourceFinder, make_2dgaussian_kernel
from photutils.utils.exceptions import NoDetectionsWarning

PHOTUTILS_GE_3 = minversion(photutils, "2.3.1.dev")

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
    convolved_data, snr_threshold, n_pixels, bkg_rms, deblend=False, mask=None
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

    Returns
    -------
    segment_img : `SegmentationImage`
        The segmentation image.
    """
    with warnings.catch_warnings():
        # suppress NoDetectionsWarning from photutils
        warnings.filterwarnings("ignore", category=NoDetectionsWarning)

        finder = SourceFinder(n_pixels, deblend=deblend, contrast=1e-4)
        threshold = snr_threshold * bkg_rms
        segment_img = finder(convolved_data, threshold, mask=mask)

        if segment_img is None:
            n_sources = 0
        else:
            if PHOTUTILS_GE_3:
                n_sources = segment_img.n_labels
            else:
                n_sources = segment_img.nlabels
        log.info(f"Detected {n_sources} sources")

    return segment_img
