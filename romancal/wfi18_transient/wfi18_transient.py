import logging
import warnings

import numpy as np
from astropy.stats import sigma_clipped_stats
from roman_datamodels import dqflags
from scipy import optimize

__all__ = ["correct_anomaly", "mask_affected_rows"]

log = logging.getLogger(__name__)


def mask_affected_rows(groupdq):
    """
    Mask the rows most affected by the transient anomaly.

    Rows masked with DO_NOT_USE are 0-1000, at the bottom of the array.

    Parameters
    ----------
    groupdq : `~numpy.ndarray`
        The group DQ image.  Updated in place.
    """
    groupdq[0, :1000, :] |= dqflags.group.DO_NOT_USE


def _frame_read_times(sca, frame_time):
    """
    Compute the read times for a single frame.

    This is a placeholder function that assumes a uniform read
    across each channel within the frame time.  A more careful
    treatment will need to account for time spent reading out the
    guide window.

    Data shape for the frame is assumed to be 4096 x 4096, with
    32 channels along the columns.

    Parameters
    ----------
    sca : int
        Detector number.
    frame_time : float
        The frame time for the exposure, in seconds.

    Returns
    -------
    read_times : `~numpy.ndarray`
        A 4096 x 4096 array containing the read time in seconds for each pixel.
    """
    nchannel = 32
    nrow = 4096
    ncol = 128
    channel_read_times = (
        np.arange(nrow * ncol).reshape(nrow, ncol) / (nrow * ncol) * frame_time
    )
    read_times = np.tile(channel_read_times, (1, nchannel))

    # Apply science -> detector flipping for read order
    if (sca % 3) == 0:
        read_times = read_times[:, ::-1]
    else:  # pragma: no cover
        # This is not expected for WFI18, but included here for completeness.
        read_times = read_times[::-1, :]
    return read_times


def _resultant_read_times(read_pattern, frame_time):
    """
    Compute the mean read time for each resultant in the exposure.

    Parameters
    ----------
    read_pattern : list of list of int
        The read pattern for the exposure, from `meta.exposure.read_pattern`.
    frame_time : float
        The frame time for the exposure, in seconds.

    Returns
    -------
    `~numpy.ndarray`
        One float value for each resultant, containing the mean read time
        in seconds.
    """
    t_resultant = []
    for pattern in read_pattern:
        group_read = [r * frame_time for r in pattern]
        t_resultant.append(np.mean(group_read))
    return np.array(t_resultant)


def _double_exp(t, a, b, tau_a, tau_b):
    """
    Evaluate a double exponential function.

    Functional form is ``a * exp(-t/tau_a) + b * exp(-t/tau_b)``.

    Parameters
    ----------
    t : float
        Time (independent variable).
    a : float
        Amplitude for the first exponential.
    b : float
        Amplitude for the second exponential.
    tau_a : float
        Time constant for the first exponential.
    tau_b : float
        Time constant for the second exponential.

    Returns
    -------
    float
        The double exponential function evaluated at time ``t``.
    """
    return a * np.exp(-t / tau_a) + b * np.exp(-t / tau_b)


def correct_anomaly(input_model, mask_rows=False):
    """
    Correct the transient first read anomaly.

    If the input model has less than 5 resultants or the fitting
    process fails, the most affected rows in the first read will be
    masked instead of fitting and removing the anomaly.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        Input model to correct. Updated in place.
    mask_rows : bool, optional
        If True, use the fallback option of just masking the
        most affected rows instead of fitting and removing
        the anomaly.

    Returns
    -------
    input_model : `~roman_datamodels.datamodels.RampModel`
        The updated model.
    """
    # Fallback option: just mask the affected rows in the first read
    if mask_rows or input_model.data.shape[0] < 5:
        log.info("Masking affected rows")
        mask_affected_rows(input_model.groupdq)
        return input_model

    # Get the read times for all pixels and resultants
    sca = input_model.meta.exposure.sca_number
    frame_time = input_model.meta.exposure.frame_time
    read_pattern = input_model.meta.exposure.read_pattern
    t_pixel = _frame_read_times(sca, frame_time)
    t_resultant = _resultant_read_times(read_pattern, frame_time)

    # Input data without reference pixels
    data = input_model.data[:, 4:-4, 4:-4]
    t_pixel = t_pixel[4:-4, 4:-4]

    # Average read time per row
    t_row = np.mean(t_pixel, axis=1)

    # The anomalous recorded value for the first read
    first_read = data[0].copy()

    # Use the 2nd and 3rd resultants to estimate the true value for the first read,
    # extrapolating backward from the count rate.
    diff_time = np.diff(t_resultant)
    estimated_rate = (data[2] - data[1]) / diff_time[1]
    estimated_first_read = data[1] - estimated_rate * diff_time[0]

    # Residual from the estimated value
    residual = first_read - estimated_first_read

    # Identify weakly illuminated pixels in each row for fitting.
    # Use a different pair of reads, to minimize bias from covariances.
    diff_data = data[4] - data[3]
    high_illum_mask = diff_data > np.median(diff_data, axis=1)

    # Compute the sigma clipped mean of the residual for each row
    # from the low illumination pixels
    residual_rows, _, _ = sigma_clipped_stats(
        residual, mask=high_illum_mask, sigma=3.0, axis=1
    )

    # Fit the sum of two exponential functions to the residual rows
    # Starting guess is here is empirical and may need adjustment.
    p0 = [residual_rows[0] * 1.1, residual_rows[0] * -0.1, t_row[300], t_row[1000]]

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                "Covariance of the parameters could not be estimated",
                optimize.OptimizeWarning,
            )
            result = optimize.curve_fit(_double_exp, t_row, residual_rows, p0=p0)
    except (ValueError, RuntimeError):
        # Fit failed - fall back on masking the data
        log.warning("Transient fit failed; masking affected rows instead.")
        mask_affected_rows(input_model.groupdq)
        return input_model

    # Calculate the correction by evaluating the model at the pixel readout times.
    a, b, tau_a, tau_b = result[0]
    log.debug(f"Fit parameters: a={a}, b={b}, tau_a={tau_a}, tau_b={tau_b}")
    correction = _double_exp(t_pixel, a, b, tau_a, tau_b)

    # Subtract the correction in the array view
    data[0] -= correction

    return input_model
