import logging
import warnings

import numpy as np
from astropy.convolution import convolve_fft
from astropy.utils.exceptions import AstropyUserWarning
from photutils.segmentation import make_2dgaussian_kernel

from romancal.datamodels import ModelLibrary

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_gaussian_kernel(fwhm):
    """
    Make a normalized 2D circular Gaussian kernel.

    The kernel will have odd sizes in both X and Y, be centered in the
    central pixel, and normalized to sum to 1.

    Parameters
    ----------
    fwhm : float
        The full-width at half-maximum (FWHM) in pixels of the 2D
        Gaussian kernel.

    Returns
    -------
    kernel : `astropy.convolution.Kernel2D`
        The output Gaussian kernel.
    """
    size = np.ceil(3 * fwhm).astype(int)
    size = size + 1 if size % 2 == 0 else size  # make size be odd
    kernel = make_2dgaussian_kernel(fwhm, size=size)  # sums to 1
    return kernel.array


def make_det_image(library, kernel_fwhm):
    """
    Make a detection image from a library of models.

    The detection image is the weighted sum of the input models, where
    the weights are SED weight divided by the read noise variance.

    Parameters
    ----------
    library : `romancal.datamodels.ModelLibrary`
        The input library of models.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) in pixels of the 2D
        Gaussian kernel used to smooth the detection image.

    Returns
    -------
    detection_data : 2D `numpy.ndarray`
        The detection image data.
    """
    if not isinstance(library, ModelLibrary):
        raise TypeError("library input must be a ModelLibrary object")

    # TODO: extend to different kernel shapes beyond Gaussian?
    # would require different/additional kernel parameters
    kernel = make_gaussian_kernel(kernel_fwhm)

    log.info(f"Making detection image with kernel FWHM={kernel_fwhm}")

    with library:
        detection_data = 0.0
        detection_var = 0.0
        wht_sum = 0.0
        all_nan_mask = False
        for i, model in enumerate(library):
            # TODO: SED weights to be defined in the asn file for each
            # input filter image
            try:
                sed_weight = library.asn["products"][0]["members"][i]["sed_weight"]
            except KeyError:
                sed_weight = 1.0

            log.info(
                f"Processing model {model.meta.filename}: "
                f"filter={model.meta.instrument.optical_element}, {sed_weight=}"
            )

            # Ideally the weights should be the inverse variance of
            # all sources of noise except the source Poisson noise. The
            # closest approximation we have is var_rnoise (read noise
            # variance) -- var_rnoise is also used as the weighting in the
            # resample step. Note that model.err is the total error, which
            # includes source Poisson noise.
            wht = sed_weight / model.var_rnoise  # inverse variance
            wht_sum += np.nan_to_num(wht, copy=False, nan=0.0)

            coverage_mask = np.isnan(model.err)
            with warnings.catch_warnings():
                # Suppress warnings about any NaNs in the data
                warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
                data_conv = convolve_fft(
                    wht * model.data,
                    kernel,
                    mask=coverage_mask,
                    preserve_nan=True,
                    normalize_kernel=False,
                )
                var_conv = convolve_fft(
                    wht**2 * model.var_rnoise,
                    kernel**2,
                    mask=coverage_mask,
                    preserve_nan=True,
                    normalize_kernel=False,
                )

            detection_data += np.nan_to_num(data_conv, copy=False, nan=0.0)
            detection_var += np.nan_to_num(var_conv, copy=False, nan=0.0)
            if i == 0:
                all_nan_mask = coverage_mask
            else:
                np.logical_and(all_nan_mask, coverage_mask, out=all_nan_mask)

            library.shelve(model, modify=False)

    # pixels that are NaN in all models are set to NaN in the output
    detection_data[all_nan_mask] = np.nan
    detection_var[all_nan_mask] = np.nan

    detection_data /= wht_sum
    detection_error = np.sqrt(detection_var) / wht_sum  # std dev error

    return detection_data / (detection_error + (detection_error == 1))


def make_detection_image(library, kernel_fwhms):
    """
    Make a detection image from a library of models.

    The detection image is maximum of the detection images from each of
    the kernels.

    Parameters
    ----------
    library : `romancal.datamodels.ModelLibrary`
        The input library of models.

    kernel_fwhms : list of float
        The full-width at half-maximum (FWHM) in pixels of the 2D
        Gaussian kernels used to smooth the detection image.

    Returns
    -------
    detection_data : 2D `numpy.ndarray`
        The detection image data.

    detection_error : 2D `numpy.ndarray`
        The detection image (standard deviation) error.
    """
    log.info("Making detection image")
    det_img = -np.inf
    for kernel_fwhm in kernel_fwhms:
        img = make_det_image(library, kernel_fwhm)
        det_img = np.fmax(det_img, img)

    return det_img
