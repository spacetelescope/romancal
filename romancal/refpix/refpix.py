from __future__ import annotations

from itertools import islice
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

import numpy as np
from scipy import fft

from .data import Aligned, ChannelFFT, Coefficients, Standard, Width

# TODO:
# 5.  Implement the equivalent of the fft_interp method as a function which takes
#     in a RefPixData object and returns a RefPixData object.
# 6.  Implement the actual correction steps:
#     a. A forward fft method on the reference pixel channels.
#     b. A the application of the correction data in Fourier space.
#     c. An inverse fft method to produce the correction data in real space.
#     d. Add method to apply the correction to a RefPixData object.
# 7.  Implement the cleanup/restore to original state.
# 8.  Implement correction control object and actual correction flow.
# 9.  Implement logging.
# 10. Integrate into pipeline step.
# 11. Add regression tests using Tyler's data.


def remove_offset(standard: Standard) -> Standard:
    """
    Use linear least squares to remove any linear offset from the data.
    """
    data = standard.combined_data
    frames, rows, columns = data.shape

    # Reshape data so that it is:
    #    [frame, frame_data]
    # where frame_data is the data for a single frame
    data = data.reshape((frames, rows * columns))

    # Craate an independent variable indexed by frame and centered at zero
    indep = np.arange(frames, dtype=data.dtype)
    indep = indep - np.mean(indep)

    # Compute sums needed for linear least squares
    sx = np.sum(indep)
    sxx = np.sum(indep**2)
    sy = np.sum(data, axis=0)
    sxy = np.matmul(indep, data, dtype=data.dtype)

    # Compute the offset (y-intercept) for the fit
    offset = (sy * sxx - sx * sxy) / (frames * sxx - sx**2)

    # Apply the offset to the data and reshape to the original shapes
    data = (data - offset).reshape((frames, rows, columns))
    offset = offset.reshape(rows, columns)

    return Standard.from_combined(data, offset)


def remove_linear_trends(aligned: Aligned) -> Aligned:
    """
    Remove any trends on the frame boundary.
    """
    data = aligned.combined_data
    channels, frames, rows, columns = data.shape

    # Create an independent variable indexed by frame and centered at zero
    t = np.arange(columns * rows, dtype=data.dtype).reshape((rows, columns))
    t = t - np.mean(t)

    # Locate the top and bottom reference pixel rows
    REF_ROWS = [*np.arange(Width.REF), *(rows - np.arange(Width.REF) - 1)[::-1]]

    # Restrict to the reference pixels and non-zero values
    t_ref = t[REF_ROWS, :]
    not_zero = data != 0

    # Fit perform a fit on a per-channel, per-frame basis
    for chan in range(channels):
        for frame in range(frames):
            # Find the non-zero reference pixels for this channel and frame
            mask = not_zero[chan, frame, REF_ROWS, :]

            # Compute the independent and dependent variables for the fit
            x_vals = t_ref[mask]
            y_vals = data[chan, frame, REF_ROWS, :][mask]

            # Skip if there is no data
            if x_vals.size < 1 or y_vals.size < 1:
                continue

            # Perform the fit using a 1st order polynomial
            m, b = np.polyfit(x_vals, y_vals, 1)

            # Remove the fit from the data
            data[chan, frame, :, :] -= (t * m + b) * not_zero[chan, frame, :, :]

    return Aligned.from_combined(data, aligned.offset)


def amp33_cosine_interpolation(aligned: Aligned) -> Aligned:
    """
    Perform Cosine weighted interpolation on the zero values of the amp33 channels
    """
    data = aligned.amp33
    channels, _, rows, columns = aligned.data.shape

    interp = np.sin(np.arange(1, channels + 1, dtype=data.dtype) * np.pi / channels)

    for frame in data:
        kern = frame.reshape(rows * columns)

        mask = (kern != 0).astype(np.int64)

        cov = np.convolve(kern, interp, mode="same")
        mask_conv = np.convolve(mask, interp, mode="same")

        kern = (cov / mask_conv).reshape(rows, columns)
        frame += kern * (1 - mask.reshape(rows, columns))

    # Fix NaN values
    np.nan_to_num(data, copy=False)

    return Aligned(aligned.data, data, aligned.offset)


def fft_interpolation(frame: np.ndarray, read_only: np.ndarray, apodize: np.ndarray):
    only = frame[read_only]

    while True:
        result = apodize * fft.rfft(frame, workers=1) / frame.size
        frame = fft.irfft(result * frame.size, workers=1)
        frame[read_only] = only

        yield frame


def amp33_fft_interpolation(aligned: Aligned, num: int = 3) -> Aligned:
    """
    FFT interpolate the zero values of the amp33 channels
    """
    frames, rows, columns = aligned.amp33.shape
    length = rows * columns

    mask = np.zeros((rows, columns), dtype=bool)
    mask[:, : -Width.PAD] = True
    mask = mask.flatten()

    read_only = np.where(mask)[0]

    data = aligned.amp33.reshape(frames, length)
    apodize = (
        1 + np.cos(2 * np.pi * np.abs(np.fft.rfftfreq(length, 1 / length)) / length)
    ) / 2

    for index, frame in enumerate(data):
        data[index, :] = next(
            islice(fft_interpolation(frame, read_only, apodize), num - 1, None)
        )

    return Aligned(aligned.data, data.reshape(frames, rows, columns), aligned.offset)


def correction(channel_fft: ChannelFFT, coeffs: Coefficients) -> np.ndarray:
    channels, columns = coeffs.gamma.shape

    # We need to multiply by 2 because we are only using half of the FFT because
    # the data is real
    normalization = columns * 2

    for idx in range(channels):
        correct = (
            np.multiply(channel_fft.left, coeffs.gamma[idx])
            + np.multiply(channel_fft.right, coeffs.zeta[idx])
            + np.multiply(channel_fft.amp33, coeffs.alpha[idx])
        ) * normalization

        correct = fft.irfft(correct).real
        yield correct

    # Add zeros in for the amp33 channel as it does not get changed
    yield np.zeros(correct.shape)


def apply_correction(aligned: Aligned, coeffs: Coefficients) -> Aligned:
    data = aligned.combined_data
    channel_fft = ChannelFFT.from_aligned(aligned)

    return channel_fft

    # Expand the correction iteratior
    correct = np.array(list(correction(channel_fft, coeffs)))
    return correct
    correct = correct.reshape(data.shape)

    return correct

    # Apply the actual correction (throw padding away)
    data[:, :, :, : -Width.PAD] -= correct.reshape(data.shape)[:, :, :, : -Width.PAD]

    return Aligned.from_combined(data, aligned.offset)
