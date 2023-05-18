from __future__ import annotations

import abc
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from roman_datamodels.datamodels import RampModel

from dataclasses import dataclass
from enum import IntEnum

import numpy as np
from astropy import units as u


class Width(IntEnum):
    REF = 4
    CHANNEL = 128
    PAD = 12


def _extract_value(data):
    if isinstance(data, u.Quantity):
        if data.unit != u.DN:
            raise ValueError(f"Input data must be in units of DN, not {data.unit}")
        data = data.value

    return data


_offset = np.ndarray | None


@dataclass
class Base(abc.ABC):
    data: np.ndarray
    amp33: np.ndarray
    offset: _offset = None

    @abc.abstractproperty
    def left(self) -> np.ndarray:
        """Left reference pixels."""
        ...

    @abc.abstractproperty
    def right(self) -> np.ndarray:
        """Right reference pixels."""
        ...

    @abc.abstractproperty
    def combined_data(self) -> np.ndarray:
        """All data in one array. (regular with amp33)"""
        ...

    @abc.abstractclassmethod
    def from_combined(cls, data: np.ndarray, offset: _offset = None) -> Base:
        """
        Extract the data and amp33 from the combined data.
        """
        ...


@dataclass
class Standard(Base):
    """
    Hold all the data in the normal detector arrangement.
    """

    @property
    def left(self) -> np.ndarray:
        return self.data[:, :, : Width.REF]

    @property
    def right(self) -> np.ndarray:
        return self.data[:, :, -Width.REF :]

    @property
    def combined_data(self) -> np.ndarray:
        return np.dstack([self.data, self.amp33])

    @classmethod
    def from_combined(cls, data: np.ndarray, offset: _offset = None) -> Standard:
        return cls(data[:, :, : -Width.CHANNEL], data[:, :, -Width.CHANNEL :], offset)

    @property
    def aligned_data(self) -> np.ndarray:
        """
        Return the aligned and padded data.
        """
        data = self.combined_data
        frames, rows, columns = data.shape
        channels = columns // Width.CHANNEL

        # First split data into channels
        #    [frame, row, channel, channel_column]
        # Second reorder so channels are first (better memory alignment)
        #    [channel, frame, row, channel_column]
        data = data.reshape((frames, rows, channels, Width.CHANNEL)).transpose(
            (2, 0, 1, 3)
        )

        # Reverse every other channel column so columns are in reading order
        data[1::2, :, :, :] = data[1::2, :, :, ::-1]

        # pad channels with zeros to account for the pause at the end of each read
        return np.pad(data, ((0, 0), (0, 0), (0, 0), (0, Width.PAD)), constant_values=0)

    @classmethod
    def from_aligned(cls, ref_data: Aligned) -> Standard:
        """
        Construct a Standard object from an Aligned object.
        """

        return cls.from_combined(ref_data.standard_data, ref_data.offset)

    @classmethod
    def from_datamodel(cls, datamodel: RampModel) -> Standard:
        data = _extract_value(datamodel.data)

        # Extract amp33
        amp33 = _extract_value(datamodel.amp33)
        # amp33 is normally a uint16, but the computation assumes it is a float
        amp33 = amp33.astype(data.dtype)

        return cls(data, amp33)


@dataclass
class Aligned(Base):
    """
    Hold all the data in the aligned and padded detector arrangement.
    """

    @property
    def left(self):
        return self.data[0, :, :, :]

    @property
    def right(self):
        return self.data[-1, :, :, :]

    @property
    def combined_data(self) -> np.ndarray:
        return np.append(self.data, [self.amp33], axis=0)

    @classmethod
    def from_combined(cls, data: np.ndarray, offset: _offset = None) -> Aligned:
        return cls(data[:-1, :, :, :], data[-1, :, :, :], offset)

    @property
    def standard_data(self) -> np.ndarray:
        """
        Return the standard data.
        """

        # Remove the padding
        data = self.combined_data[:, :, :, : -Width.PAD]
        channels, frames, rows, columns = data.shape

        # Reverse every other channel column so columns are in reading order
        data[1::2, :, :, :] = data[1::2, :, :, ::-1]

        # Undo channel reordering to:
        #    [frame, row, channel, channel_column]
        # Second combine the channels into columns:
        #    [frame, row, column]
        return data.transpose((1, 2, 0, 3)).reshape((frames, rows, columns * channels))

    @classmethod
    def from_standard(cls, ref_data: Standard) -> Aligned:
        """
        Construct an Aligned object from a Standard object.
        """

        return cls.from_combined(ref_data.aligned_data, ref_data.offset)
