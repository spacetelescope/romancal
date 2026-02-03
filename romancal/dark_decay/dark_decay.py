import numpy as np

from romancal.lib.basic_utils import frame_read_times

__all__ = ["subtract_dark_decay"]


def subtract_dark_decay(data, amplitude, time_constant, frame_time, read_pattern, sca):
    """
    Subtract dark decay correction in place.

    Parameters
    ----------
    data : np.ndarray
        Data array to correct. Updated in place.
    amplitude : array-like
        Decay amplitude for the detector from the reference file.
    time_constant : array-like
        Decay time constant for the detector from the reference file.
    frame_time : float
        The frame time for the exposure, in seconds.
    read_pattern : list of list of int
        The read pattern for the exposure.
    sca : int
        The number of the WFI detector (1-18).
    """
    for i in range(data.shape[0]):  # loop over resultants
        corrections = [
            amplitude * np.exp(-frame_read_times(frame_time, sca, j) / time_constant)
            for j in read_pattern[i]
        ]
        data[i] -= np.mean(corrections, axis=0)
