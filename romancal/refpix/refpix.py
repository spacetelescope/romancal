from __future__ import annotations

from itertools import islice
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from roman_datamodels import RampModel

from dataclasses import dataclass
from enum import IntEnum, StrEnum, auto

import numpy as np
from astropy import units as u
from scipy import fft

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


class Width(IntEnum):
    REF = 4
    CHANNEL = 128
    PAD = 12


class Arrangement(StrEnum):
    STANDARD = auto()
    ALIGNED = auto()


def _extract_value(data):
    if isinstance(data, u.Quantity):
        if data.unit != u.DN:
            raise ValueError(f"Input data must be in units of DN, not {data.unit}")
        data = data.value

    return data


@dataclass
class RefPixData:
    data: np.ndarray
    amp33: np.ndarray
    offset: np.ndarray = None
    arrangement: Arrangement = Arrangement.STANDARD

    @classmethod
    def from_datamodel(cls, datamodel: RampModel):
        data = _extract_value(datamodel.data)

        # Extract amp33
        amp33 = _extract_value(datamodel.amp33)
        # amp33 is normally a uint16, but the computation assumes it is a float
        amp33 = amp33.astype(data.dtype)

        return cls(data, amp33)

    @property
    def combine_data(self):
        if self.arrangement == Arrangement.STANDARD:
            return np.dstack([self.data, self.amp33])
        else:
            data = self.split_channels
            channels, frames, rows, columns = data.shape

            # Undo channel reordering to:
            #    [frame, row, channel, channel_column]
            # Second combine the channels into columns:
            #    [frame, row, column]
            return data.transpose((1, 2, 0, 3)).reshape(
                (frames, rows, columns * channels)
            )

    @property
    def left(self):
        if self.arrangement == Arrangement.STANDARD:
            return self.data[:, :, : Width.REF]
        else:
            return self.data[0, :, :, :]

    @property
    def right(self):
        if self.arrangement == Arrangement.STANDARD:
            return self.data[:, :, -Width.REF :]
        else:
            return self.data[-1, :, :, :]

    @classmethod
    def from_combined_data(cls, data, offset=None):
        return cls(data[:, :, : -Width.CHANNEL], data[:, :, -Width.CHANNEL :], offset)

    @property
    def split_channels(self):
        if self.arrangement == Arrangement.STANDARD:
            data = self.combine_data
            frames, rows, columns = data.shape
            channels = columns // Width.CHANNEL

            # First split data into channels
            #    [frame, row, channel, channel_column]
            # Second reorder so channels are first (better memory alignment)
            #    [channel, frame, row, channel_column]
            return data.reshape((frames, rows, channels, Width.CHANNEL)).transpose(
                (2, 0, 1, 3)
            )
        else:
            return np.append(self.data, [self.amp33], axis=0)

    @property
    def aligned_channels(self):
        if self.arrangement == Arrangement.ALIGNED:
            return self.split_channels
        else:
            data = self.split_channels

            # Reverse every other channel column so columns are in reading order
            data[1::2, :, :, :] = data[1::2, :, :, ::-1]

            # pad channels with zeros to account for the pause at the end of each read
            return np.pad(
                data, ((0, 0), (0, 0), (0, 0), (0, Width.PAD)), constant_values=0
            )

    @classmethod
    def from_split(cls, data, offset=None):
        return cls(data[:-1, :, :, :], data[-1, :, :, :], offset, Arrangement.ALIGNED)


def remove_offset(ref_data: RefPixData) -> RefPixData:
    """
    Use linear least squares to remove any linear offset from the data.
    """
    data = ref_data.combine_data
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

    return RefPixData.from_combined_data(data, offset)


def remove_linear_trends(ref_data: RefPixData) -> RefPixData:
    """
    Remove any trends on the frame boundary.
    """
    data = ref_data.aligned_channels
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

    return RefPixData.from_split(data, ref_data.offset)


def cosine_interpolate_amp33(ref_data: RefPixData) -> RefPixData:
    """
    Perform Cosine weighted interpolation on the zero values of the amp33 channels
    """

    if ref_data.arrangement != Arrangement.ALIGNED:
        ref_data = RefPixData.from_split(ref_data.aligned_channels, ref_data.offset)

    data = ref_data.amp33
    channels, _, rows, columns = ref_data.data.shape

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

    return RefPixData(ref_data.data, data, ref_data.offset, ref_data.arrangement)


def fft_interpolate(frame, read_only, apodize):
    data = frame[read_only]
    while True:
        result = apodize * fft.rfft(data, workers=1) / data.size
        data = fft.ifft(result * data.size, workers=1)

        yield data


def fft_interpolate_amp22(ref_data: RefPixData, num: int = 50) -> RefPixData:
    """
    FFT interpolate the zero values of the amp22 channels
    """

    frames, rows, columns = ref_data.amp33.shape
    length = rows * columns

    mask = np.zeros(columns, dtype=bool)
    mask[: -Width.PAD] = True
    read_only = np.where(mask)

    data = ref_data.amp33.reshape(frames, length)
    apodize = (
        1 + np.cos(2 * np.pi * np.abs(np.fft.rfftfreq(length, 1 / length)) / length)
    ) / 2

    for frame in data:
        frame[read_only] = next(
            islice(fft_interpolate(frame, read_only, apodize), num, None)
        )

    return data
