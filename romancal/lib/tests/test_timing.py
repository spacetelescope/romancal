import numpy as np
import pytest

from romancal.lib.basic_utils import frame_read_times


def test_timing():

    frametime_ref = 3.16 # in seconds
    readtimes_nogw = frame_read_times(
        frametime_ref,
        1,
        frame_number=0,
        stride_gw=256,
        gw_pseudotime_factor=0
    )

    readtimes_nominal = frame_read_times(
        frametime_ref,
        1,
        frame_number=0,
        stride_gw=256,
        gw_pseudotime_factor=1
    )

    readtimes_default = frame_read_times(
        frametime_ref,
        1,
        frame_number=0,
        stride_gw=256,
        gw_pseudotime_factor=1.5
    )

    # Row-to-row differences in time
    rowdiffs_nogw = np.diff(readtimes_nogw, axis=0)[:, 0]
    rowdiffs_nominal = np.diff(readtimes_nominal, axis=0)[:, 0]
    rowdiffs_default = np.diff(readtimes_default, axis=0)[:, 0]

    # Uniform gaps between rows without the guide window
    assert np.std(rowdiffs_nogw) < 1e-6

    # Large gap at Row 256 with the guide window
    assert(rowdiffs_nominal[255] > 100*rowdiffs_nominal[254])

    # Slightly less than a factor of 1.5 between the guide window gaps
    # when using a factor of 1.5 for gw_pseudotime_factor.  This is not
    # exactly 1.5 because of the normalization to the total frame time.
    assert(np.isclose(rowdiffs_default[255]/rowdiffs_nominal[255], 1.45, 0.05))

    # Smaller time spread within rows when allocating some of the time
    # to the guide window
    assert(np.all(np.std(readtimes_nominal, axis=1) <
                  np.std(readtimes_nogw, axis=1)))

    # All frame times check out: pixel read times go from zero to
    # the reference frame time.
    for arr in [readtimes_nogw, readtimes_default, readtimes_nominal]:
        assert(np.amax(arr) <= frametime_ref)
        assert(np.amax(arr) > 0.95*frametime_ref)
        assert(np.amin(arr) < 0.05*frametime_ref)
        assert(np.amin(arr) >= 0)
        
