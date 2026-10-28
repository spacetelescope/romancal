from __future__ import annotations

import logging
from typing import TYPE_CHECKING
import numpy as np

from astropy.stats import mad_std
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
    smooth_sigma: float = 8.0,
) -> np.ndarray:
    """
    Dynamically computes row-by-row weights [Left, Right, Default, Amp33]
    using uniform-weight linear ramp fitting, Gaussian slope smoothing, 
    and Tim's continuous [0, 1] sigma-tapered validity mapping.
    """
    n_res, n_rows, _ = datamodel.data.shape

    # 1. Isolate and collapse Left and Right reference columns (shape: n_res, n_rows)
    left_refs = np.nanmedian(datamodel.data[:, :, :4], axis=2).astype(np.float32)
    right_refs = np.nanmedian(datamodel.data[:, :, -4:], axis=2).astype(np.float32)

    # 2. Vectorized Uniform Weight Linear Ramp Fit (Slopes in DN / resultant)
    x = np.arange(n_res, dtype=np.float32)
    x_diff = x - np.mean(x)
    x_ssq = np.sum(x_diff ** 2)

    left_temp_diff = left_refs - np.nanmean(left_refs, axis=0)[np.newaxis, :]
    left_slopes = np.sum(x_diff[:, np.newaxis] * left_temp_diff, axis=0) / x_ssq

    right_temp_diff = right_refs - np.nanmean(right_refs, axis=0)[np.newaxis, :]
    right_slopes = np.sum(x_diff[:, np.newaxis] * right_temp_diff, axis=0) / x_ssq

    # 3. Spatial Median Centering across the 4096 rows
    left_dev = left_slopes - np.nanmedian(left_slopes)
    right_dev = right_slopes - np.nanmedian(right_slopes)

    # 4A. Tim's signed slope smoothing (catches coherent one-sided drifts)
    left_smooth = gaussian_filter1d(left_dev, sigma=smooth_sigma)
    right_smooth = gaussian_filter1d(right_dev, sigma=smooth_sigma)

    left_noise = max(1.5 * float(np.nanmedian(np.abs(left_smooth))), 1e-12)
    right_noise = max(1.5 * float(np.nanmedian(np.abs(right_smooth))), 1e-12)

    left_snr_signed = np.abs(left_smooth) / left_noise
    right_snr_signed = np.abs(right_smooth) / right_noise

    # 4B. Rectified slope envelope (prevents bipolar yellow/blue stripes in the wings
    # and core from canceling out during smoothing, normalized by its own MAD noise)
    left_env = gaussian_filter1d(np.abs(left_dev), sigma=smooth_sigma)
    right_env = gaussian_filter1d(np.abs(right_dev), sigma=smooth_sigma)

    left_env_med = float(np.nanmedian(left_env))
    right_env_med = float(np.nanmedian(right_env))

    left_env_noise = max(1.5 * float(np.nanmedian(np.abs(left_env - left_env_med))), 1e-12)
    right_env_noise = max(1.5 * float(np.nanmedian(np.abs(right_env - right_env_med))), 1e-12)

    left_snr_env = (left_env - left_env_med) / left_env_noise
    right_snr_env = (right_env - right_env_med) / right_env_noise

    # Combine signed and rectified SNR so neither core sign-flips nor wings drop out
    left_snr = np.maximum(left_snr_signed, left_snr_env)
    right_snr = np.maximum(right_snr_signed, right_snr_env)

    # 5. Tim's continuous validity weight: tapers from 1 (valid) to 0 (corrupted)
    # across [sigma_thresh - 0.5, sigma_thresh + 0.5]
    w_left_valid = np.clip((sigma_thresh + 0.5) - left_snr, 0.0, 1.0)
    w_right_valid = np.clip((sigma_thresh + 0.5) - right_snr, 0.0, 1.0)

    c_left = 1.0 - w_left_valid
    c_right = 1.0 - w_right_valid

    log.info(
        f"Dynamic weights active on {np.sum(c_left > 0)} Left rows "
        f"(full weight on {np.sum(c_left == 1.0)} rows)."
    )
    log.info(
        f"Dynamic weights active on {np.sum(c_right > 0)} Right rows "
        f"(full weight on {np.sum(c_right == 1.0)} rows)."
    )

    # 6. Bilinearly map valid/corrupted weights to [Left, Right, Default, Amp33]
    weights = np.zeros((n_rows, 4), dtype=np.float32)
    weights[:, 0] = w_left_valid * c_right          # Left-only (Left valid, Right corrupted)
    weights[:, 1] = c_left * w_right_valid          # Right-only (Left corrupted, Right valid)
    weights[:, 2] = w_left_valid * w_right_valid    # Default (Both valid)
    weights[:, 3] = c_left * c_right                # Amp33 (Both corrupted)

    np.save(
        '/grp/roman/hgk/romcom/crosstalk_esque_issue/dr15_dynamic_weights_final.npy',
        weights,
    )

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
    # Substituting R <- L requires gamma_left = gamma + zeta so healthy baseline drift
    # is still fully subtracted when one side is zeroed out.
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
        weights = _compute_dynamic_weights(datamodel, sigma_thresh=sigma_thresh, smooth_sigma=smooth_sigma)
        w_left = weights[:, 0][:, np.newaxis]
        w_right = weights[:, 1][:, np.newaxis]
        w_default = weights[:, 2][:, np.newaxis]
        w_amp33 = weights[:, 3][:, np.newaxis]

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

            standard.data[:] = (
                (left_data * w_left)
                + (right_data * w_right)
                + (default_data * w_default)
                + (amp33_data * w_amp33)
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