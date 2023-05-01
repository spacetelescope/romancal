from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from roman_datamodels import RampModel

from dataclasses import dataclass
from enum import IntEnum, StrEnum, auto

import numpy as np
from astropy import units as u

# TODO:
# 1.  Add offset/reversed channel tracking and undo to the RefPixData object.
# 3.  Refactor the remove_linear_trends method into a function which takes in
#     a RefPixData object and returns a RefPixData object.
# 4.  Implement the equivalent of the interp_zeros_channel_fun as a function which
#     takes in a RefPixData object and returns a RefPixData object.
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
    SPLIT = auto()


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
        data = self.split_channels

        # Reverse every other channel column so columns are in reading order
        data[1::2, :, :, :] = data[1::2, :, :, ::-1]

        # pad channels with zeros to account for the pause at the end of each read
        return np.pad(data, ((0, 0), (0, 0), (0, 0), (0, Width.PAD)), constant_values=0)

    @classmethod
    def from_split(cls, data, offset=None):
        return cls(data[:-1, :, :, :], data[-1, :, :, :], offset, Arrangement.SPLIT)

    def remove_linear_trends(self):
        data = self.aligned_channels
        channels, frames, rows, columns = data.shape

        t = np.arange(columns * rows, dtype=data.dtype).reshape((rows, columns))
        t = t - np.mean(t)

        REF_ROWS = [*np.arange(Width.REF), *(rows - np.arange(Width.REF) - 1)[::-1]]

        t_ref = t[REF_ROWS, :]
        not_zero = data != 0

        for chan in range(channels):
            for frame in range(frames):
                mask = not_zero[chan, frame, REF_ROWS, :]

                x_vals = t_ref[mask]
                y_vals = data[chan, frame, REF_ROWS, :][mask]

                if x_vals.size < 1 or y_vals.size < 1:
                    continue

                m, b = np.polyfit(x_vals, y_vals, 1)
                data[chan, frame, :, :] -= (t * m + b) * not_zero[chan, frame, :, :]

        return data


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
