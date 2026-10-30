from __future__ import annotations

import logging
from typing import TYPE_CHECKING
import numpy as np

from scipy.ndimage import gaussian_filter1d
from roman_datamodels.dqflags import pixel
from romancal.ramp_fitting.ramp_fit_step import slopes_uniform_weights

if TYPE_CHECKING:
    from roman_datamodels.datamodels import RampModel, RefpixRefModel

from ._data import Coefficients, StandardView

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _compute_dynamic_weights(
    datamodel: RampModel,
    sigma_thresh: float = 3.5,
    sigma_width: float = 1.0,
    smooth_sigma: float = 12.0,
) -> dict:
    """
    We have four different kinds of fits---using all of the border reference pixels, 
    using only the left reference pixels, using only the right reference pixels, 
    and using no border reference pixels. We want to weight these 4 kinds of fits 
    according to how consistent the border pixels are with typical fluctuations in the border. 
    We derive left and right weights between 0 and 1, and then get 4 weights that sum to one 
    and correspond to bilinear interpolation among these four options along the left and 
    right weight axes.
    """

    # 1. Isolate and collapse Left and Right reference columns
    # Using nanmean to match fundamental refpix correction math
    left_refs = np.nanmean(datamodel.data[:, :, :4], axis=2).astype(np.float32)
    right_refs = np.nanmean(datamodel.data[:, :, -4:], axis=2).astype(np.float32)

    # 2. Vectorized Uniform Weight Linear Ramp Fit (Slopes in DN / resultant)
    # This is the maximum likelihood slope in the read noise limit.
    ml_weights = slopes_uniform_weights(datamodel, return_weights=True)

    left_slopes = np.sum(left_refs * ml_weights[:, None], axis=0)
    right_slopes = np.sum(right_refs * ml_weights[:, None], axis=0)

    # 3. Signed slope smoothing (catches coherent one-sided drifts)
    # Spatial median centering is mathematically required. Healthy ref pixels 
    # accumulate a positive baseline drift over the ramp. Without this, the global 
    # DC drift inflates the MAD noise floor and results in 0 detections.
    left_dev = left_slopes - np.nanmedian(left_slopes)
    right_dev = right_slopes - np.nanmedian(right_slopes)

    left_smooth = gaussian_filter1d(left_dev, sigma=smooth_sigma)
    right_smooth = gaussian_filter1d(right_dev, sigma=smooth_sigma)

    # The factor of 1.5 is the (approximate) ratio of the median absolute deviation
    # to the standard deviation assuming Gaussianity.
    # The 1e-12 floor is intended to do no harm in the all zero / constant / 
    # reference-pixel-not-simulated case.
    left_noise = max(1.5 * np.nanmedian(np.abs(left_smooth)), 1e-12)
    right_noise = max(1.5 * np.nanmedian(np.abs(right_smooth)), 1e-12)

    left_snr = np.abs(left_smooth) / left_noise
    right_snr = np.abs(right_smooth) / right_noise

    # 4. Continuous validity weight: tapers from 1 (valid) to 0 (corrupted)
    # across [sigma_thresh - (sigma_width/2), sigma_thresh + (sigma_width/2)]
    # For example, with default sigma_thresh = 3.5, sigma_width = 1.0, 
    # 4-sigma outliers get no weight while 3-sigma outliers get full weight.
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
    weight_mode: str = 'dynamic',
    sigma_thresh: float = 3.5,
    sigma_width: float = 1.0,
    smooth_sigma: float = 12.0,
) -> RampModel:
    """
    Organize the steps to run the reference pixel correction.
    Apart from the offset, the correction is independent between resultants, so
    the ramp is corrected one resultant at a time. The intermediate channel view
    and correction arrays are each larger than the ramp itself, so doing this
    keeps the peak memory usage proportional to a single resultant rather than
    to the whole ramp. The offset is fit over all of the resultants, so it has to be accumulated
    # before any of them are corrected.
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
    if weight_mode == 'dynamic':
        weights = _compute_dynamic_weights(
            datamodel,
            sigma_thresh=sigma_thresh,
            sigma_width=sigma_width,
            smooth_sigma=smooth_sigma
        )
    elif weight_mode != 'static':
        # Scaffolding for unit tests forcing specific states ('left', 'right', 'default', 'none')
        n_rows = datamodel.data.shape[1]
        weights = {
            'left_only': np.zeros(n_rows, dtype=np.float32),
            'right_only': np.zeros(n_rows, dtype=np.float32),
            'all_valid': np.zeros(n_rows, dtype=np.float32),
            'amp33_only': np.zeros(n_rows, dtype=np.float32)
        }
        if weight_mode == 'left':
            weights['left_only'] += 1.0
        elif weight_mode == 'right':
            weights['right_only'] += 1.0
        elif weight_mode == 'default':
            weights['all_valid'] += 1.0
        elif weight_mode == 'none':
            weights['amp33_only'] += 1.0
        else:
            raise ValueError(f"Unsupported weight_mode: {weight_mode}")
            

    offset = None
    if remove_offset:
        offset = StandardView.compute_offset(datamodel.data, datamodel.amp33)
        log.debug("Computed the general offset of the data, to be re-applied later.")

    log.debug(
        "Correcting the data one resultant at a time "
        f"(remove_trends={remove_trends}, cosine_interpolate={cosine_interpolate}, "
        f"fft_interpolate={fft_interpolate})."
    )

    def _get_corrected_state(standard_view, coeffs_state, orig_data):
        """Helper to apply a specific coefficient state to a single resultant."""
        standard_view.data[:] = orig_data.copy()
        ch = standard_view.channels
        if remove_trends:
            ch = ch.remove_trends()
        if cosine_interpolate:
            ch = ch.cosine_interpolate()
        if fft_interpolate:
            ch = ch.fft_interpolate()
        
        # Perform the reference pixel correction
        return ch.apply_correction(coeffs_state).data.copy()

    # The ramp is corrected one resultant at a time. Reusing _get_corrected_state() 
    # for each coefficient state ensures peak memory usage remains proportional to a 
    # single resultant rather than the entire ramp.
    for resultant in range(datamodel.data.shape[0]):
        # Read a single resultant from the datamodel into the standard view
        standard = StandardView.from_datamodel(datamodel, resultant)

        if offset is not None:
            # Remove offset from the data
            standard.remove_offset(offset)

        orig_data = standard.data.copy()

        if weight_mode == 'static':
            # Run the single default correction pipeline without blending overhead
            standard.data[:] = _get_corrected_state(standard, coeffs_default, orig_data)
        else:
            # Generate the 4 coefficient states and blend via dynamic validity weights
            left_data = _get_corrected_state(standard, coeffs_left, orig_data)
            right_data = _get_corrected_state(standard, coeffs_right, orig_data)
            amp33_data = _get_corrected_state(standard, coeffs_amp33, orig_data)
            default_data = _get_corrected_state(standard, coeffs_default, orig_data)

            # Dictionary unpacked directly into the blend
            standard.data[:] = (
                (left_data * weights['left_only'][:, None])
                + (right_data * weights['right_only'][:, None])
                + (default_data * weights['all_valid'][:, None])
                + (amp33_data * weights['amp33_only'][:, None])
            )

        if offset is not None:
            standard.apply_offset()
            
        # Write the resultant back to the datamodel

        standard.update(datamodel, resultant)

    log.debug("Updated the datamodel with the corrected data.")

    return datamodel


def _mask_bad_ref_pixels(datamodel):
    """
    Zero out BAD_REF_PIXEL pixels before refpix correction.
    """
    dq = datamodel.pixeldq
    bad_ref_mask = (dq & pixel.BAD_REF_PIXEL) != 0

    # Protect against incorrectly shaped mock data in the unit tests.
    # Production datamodels enforce strict matching schemas.
    if bad_ref_mask.shape == datamodel.data.shape[-len(bad_ref_mask.shape) :]:
        datamodel.data[..., bad_ref_mask] = 0
    else:
        log.debug(
            "Skipping bad_ref_pixel masking: pixeldq shape %s does not match data trailing shape",
            bad_ref_mask.shape,
        )

    return datamodel