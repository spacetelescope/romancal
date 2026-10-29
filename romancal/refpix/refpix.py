from __future__ import annotations

import logging
from typing import TYPE_CHECKING
import numpy as np

from scipy.ndimage import gaussian_filter1d
from roman_datamodels.dqflags import pixel

if TYPE_CHECKING:
    from roman_datamodels.datamodels import RampModel, RefpixRefModel

from ._data import Coefficients, StandardView

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _compute_dynamic_weights(
    datamodel: RampModel,
    sigma_thresh: float = 2.5,
    sigma_width: float = 1.0,
    smooth_sigma: float = 8.0,
) -> dict:
    """
    Dynamically computes row-by-row weights for the 4 reference pixel states.
    Uses max-likelihood read-noise slope fitting, Gaussian slope smoothing, 
    and a continuous sigma-tapered validity mapping.
    """
    n_res, n_rows, _ = datamodel.data.shape

    # 1. Isolate and collapse Left and Right reference columns
    # Using nanmean per Tim's feedback to match fundamental refpix correction math
    left_refs = np.nanmean(datamodel.data[:, :, :4], axis=2).astype(np.float32)
    right_refs = np.nanmean(datamodel.data[:, :, -4:], axis=2).astype(np.float32)

    # 2. Vectorized Uniform Weight Linear Ramp Fit (Slopes in DN / resultant)
    # This is the maximum likelihood slope in the read noise limit.
    # read_pattern is a list of lists of the read numbers that constitute each resultant.
    ni = np.array([len(t) for t in datamodel.meta.exposure.read_pattern])
    ti = np.array([np.mean(t) for t in datamodel.meta.exposure.read_pattern])
    N = np.sum(ni)
    Nt = np.sum(ni * ti)
    Ntt = np.sum(ni * ti**2)
    
    # Calculate exact optimal weights for the slope fit
    ml_weights = (N * ni * ti - Nt * ni) / (N * Ntt - Nt**2)

    left_slopes = np.sum(left_refs * ml_weights[:, None], axis=0)
    right_slopes = np.sum(right_refs * ml_weights[:, None], axis=0)

    # 3. Tim's signed slope smoothing (catches coherent one-sided drifts)
    # RESTORED: Spatial median centering is mathematically required. Healthy ref pixels 
    # accumulate a positive baseline drift over the ramp. Without this, the global 
    # DC drift inflates the MAD noise floor and results in 0 detections.
    left_dev = left_slopes - np.nanmedian(left_slopes)
    right_dev = right_slopes - np.nanmedian(right_slopes)

    left_smooth = gaussian_filter1d(left_dev, sigma=smooth_sigma)
    right_smooth = gaussian_filter1d(right_dev, sigma=smooth_sigma)

    # The factor of 1.5 is the (approximate) ratio of the median absolute deviation
    # to the standard deviation assuming Gaussianity.
    left_noise = max(1.5 * np.nanmedian(np.abs(left_smooth)), 1e-12)
    right_noise = max(1.5 * np.nanmedian(np.abs(right_smooth)), 1e-12)

    left_snr = np.abs(left_smooth) / left_noise
    right_snr = np.abs(right_smooth) / right_noise

    # 4. Continuous validity weight: tapers from 1 (valid) to 0 (corrupted)
    # across [sigma_thresh - (sigma_width/2), sigma_thresh + (sigma_width/2)]
    w_left_valid = np.clip((sigma_thresh - left_snr) / sigma_width + 0.5, 0.0, 1.0)
    w_right_valid = np.clip((sigma_thresh - right_snr) / sigma_width + 0.5, 0.0, 1.0)

    log.info(
        f"Dynamic weights active on {np.sum((1.0 - w_left_valid) > 0)} Left rows "
        f"(full weight on {np.sum(w_left_valid == 0.0)} rows)."
    )
    log.info(
        f"Dynamic weights active on {np.sum((1.0 - w_right_valid) > 0)} Right rows "
        f"(full weight on {np.sum(w_right_valid == 0.0)} rows)."
    )

    # 5. Map valid/corrupted weights to a dictionary of the 4 blend states
    weights = {
        'left_only': w_left_valid * (1 - w_right_valid),    # left valid, right bad
        'right_only': w_right_valid * (1 - w_left_valid),   # left bad, right valid
        'all_valid': w_left_valid * w_right_valid,          # Default (both valid)
        'amp33_only': (1 - w_left_valid) * (1 - w_right_valid) # left and right both bad
    }

    return weights


def run_steps(
    datamodel: RampModel,
    refs: RefpixRefModel,
    zero_bad_ref_pix: bool = True,
    remove_offset: bool = True,
    remove_trends: bool = True,
    cosine_interpolate: bool = True,
    fft_interpolate: bool = True,
    dynamic_weights: str = 'dynamic',
    sigma_thresh: float = 2.5,
    sigma_width: float = 1.0,
    smooth_sigma: float = 8.0,
) -> RampModel:
    """
    Organize the steps to run the reference pixel correction.
    """
    if zero_bad_ref_pix:
        _mask_bad_ref_pixels(datamodel)

    # 1. Load standard CRDS coefficients (State: Default)
    coeffs_default = Coefficients.from_ref(refs)

    # 2. Generate the 3 isolated coefficient states while CONSERVING (gamma + zeta)
    coeffs_left = Coefficients(
        alpha=coeffs_default.alpha,
        gamma=coeffs_default.gamma + coeffs_default.zeta,
        zeta=np.zeros_like(coeffs_default.zeta),
    )

    coeffs_right = Coefficients(
        alpha=coeffs_default.alpha,
        gamma=np.zeros_like(coeffs_default.gamma),
        zeta=coeffs_default.zeta + coeffs_default.gamma,
    )

    coeffs_amp33 = Coefficients(
        alpha=coeffs_default.alpha,
        gamma=np.zeros_like(coeffs_default.gamma),
        zeta=np.zeros_like(coeffs_default.zeta),
    )

    # 3. Compute dynamic weights ONCE on the unmutated ramp before the resultant loop
    if dynamic_weights != 'default':
        weights = _compute_dynamic_weights(
            datamodel, 
            sigma_thresh=sigma_thresh, 
            sigma_width=sigma_width,
            smooth_sigma=smooth_sigma
        )

    offset = None
    if remove_offset:
        offset = StandardView.compute_offset(datamodel.data, datamodel.amp33)
        log.debug("Computed the general offset of the data, to be re-applied later.")

    log.debug(
        "Correcting the data one resultant at a time "
        f"(remove_trends={remove_trends}, cosine_interpolate={cosine_interpolate}, "
        f"fft_interpolate={fft_interpolate})."
    )

    for resultant in range(datamodel.data.shape[0]):
        standard = StandardView.from_datamodel(datamodel, resultant)

        if offset is not None:
            standard.remove_offset(offset)

        if dynamic_weights == 'default':
            channel = standard.channels
            if remove_trends:
                channel = channel.remove_trends()
            if cosine_interpolate:
                channel = channel.cosine_interpolate()
            if fft_interpolate:
                channel = channel.fft_interpolate()
            standard = channel.apply_correction(coeffs_default)
        else:
            original_data = standard.data.copy()

            def _get_corrected_state(coeffs_state):
                standard.data[:] = original_data.copy()
                ch = standard.channels
                if remove_trends:
                    ch = ch.remove_trends()
                if cosine_interpolate:
                    ch = ch.cosine_interpolate()
                if fft_interpolate:
                    ch = ch.fft_interpolate()
                return ch.apply_correction(coeffs_state).data.copy()

            left_data = _get_corrected_state(coeffs_left)
            right_data = _get_corrected_state(coeffs_right)
            amp33_data = _get_corrected_state(coeffs_amp33)
            default_data = _get_corrected_state(coeffs_default)

            # Dictionary unpacked directly into the blend per Tim's request
            standard.data[:] = (
                (left_data * weights['left_only'][:, None])
                + (right_data * weights['right_only'][:, None])
                + (default_data * weights['all_valid'][:, None])
                + (amp33_data * weights['amp33_only'][:, None])
            )

        if offset is not None:
            standard.apply_offset()

        standard.update(datamodel, resultant)

    log.debug("Updated the datamodel with the corrected data.")

    return datamodel


def _mask_bad_ref_pixels(datamodel):
    """
    Zero out BAD_REF_PIXEL pixels before refpix correction.
    """
    dq = datamodel.pixeldq
    bad_ref_mask = (dq & pixel.BAD_REF_PIXEL) != 0

    if bad_ref_mask.shape == datamodel.data.shape[-len(bad_ref_mask.shape) :]:
        datamodel.data[..., bad_ref_mask] = 0
    else:
        log.debug(
            "Skipping bad_ref_pixel masking: pixeldq shape %s does not match data trailing shape",
            bad_ref_mask.shape,
        )

    return datamodel