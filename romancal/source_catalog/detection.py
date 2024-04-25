"""
Module to detect sources using image segmentation.
"""

import logging
import math
import warnings

from astropy.convolution import convolve
from astropy.units import Quantity
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
    data : 2D `astropy.units.Quantity`
        The input 2D Quantity array. The data is assumed to be
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
    if not isinstance(data, Quantity):
        raise ValueError("Input model must be a Quantity array.")

    size = math.ceil(kernel_fwhm * 3)
    size = size + 1 if size % 2 == 0 else size  # make size be odd
    kernel = make_2dgaussian_kernel(kernel_fwhm, size=size)  # normalized to 1

    with warnings.catch_warnings():
        # suppress warnings caused by large masked areas
        warnings.simplefilter("ignore", AstropyUserWarning)
        return convolve(data, kernel, mask=mask)


def make_segmentation_image(
    convolved_data, snr_threshold, npixels, bkg_rms, deblend=False, mask=None
):
    """
    Make a segmentation image from a model image.

    Parameters
    ----------
    convolved_data : 2D `numpy.ndarray`
        The background-subtracted convolved data array.
    snr_threshold : float
        The per-pixel signal-to-noise ratio threshold for detection.
    npixels : int
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

        finder = SourceFinder(npixels, deblend=deblend)
        threshold = snr_threshold * bkg_rms
        segment_img = finder(convolved_data, threshold, mask=mask)

        if segment_img is None:
            nsources = 0
        else:
            nsources = segment_img.nlabels
        log.info(f"Detected {nsources} sources")

    return segment_img
