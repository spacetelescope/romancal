import numpy as np
from numpy.testing import assert_array_equal
import pytest

from jwst.datamodels import RampModel as J_RampModel

from romancal.datamodels import RampModel as R_RampModel
from romancal.datamodels import dqflags
from romancal.datamodels import GainModel, ReadNoiseModel

from ..jump import detect_jumps

# import multiprocessing # may need to enable later

from astropy.time import Time


JUMP_DET = dqflags.group["JUMP_DET"]
DO_NOT_USE = dqflags.group["DO_NOT_USE"]
GOOD = dqflags.group["GOOD"]
SATURATED = dqflags.group["SATURATED"]
NO_GAIN_VALUE = dqflags.pixel["NO_GAIN_VALUE"]


def test_nocrs_noflux(setup_inputs):
    """"
    All pixel values are zero. So slope should be zero
    """
    model, rnoise, gain = setup_inputs(ngroups=5)
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, True)

    assert np.max(out_model.groupdq) == GOOD


def test_nocrs_noflux_badgain_pixel(setup_inputs):
    """"
    all pixel values are zero. So slope should be zero, pixel with bad gain should
    have pixel dq set to 'NO_GAIN_VALUE' and 'DO_NOT_USE'
    """
    model, rnoise, gain = setup_inputs(ngroups=5, nrows=20, ncols=20)
    gain.data[7, 7] = -10  # bad gain
    gain.data[17, 17] = np.nan  # bad gain
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, True)

    # 2 bits are set for each pixel, so use bitwise_and to check is set
    assert np.bitwise_and(out_model.pixeldq[7, 7], NO_GAIN_VALUE) == NO_GAIN_VALUE
    assert np.bitwise_and(out_model.pixeldq[7, 7], DO_NOT_USE) == DO_NOT_USE
    assert np.bitwise_and(out_model.pixeldq[17, 17], NO_GAIN_VALUE) == NO_GAIN_VALUE
    assert np.bitwise_and(out_model.pixeldq[17, 17], DO_NOT_USE) == DO_NOT_USE


def test_onecr_10_groups_neighbors_flagged(setup_inputs):
    """"
    A single CR in a 10 group exposure
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10, ncols=10,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 5, 5]) == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 6] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == JUMP_DET


def test_nocr_100_groups_nframes1(setup_inputs):
    """"
    NO CR in a 100 group exposure to make sure that frames_per_group is passed correctly to
    twopoint_difference. This test recreates the problem found in issue #4571.
    """
    grouptime = 3.0
    ingain = 1
    inreadnoise = 7.0
    ngroups = 100
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nframes=1,
                                       nrows=10, ncols=10, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 14.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 27.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 38.0
    model.data[0, 5, 5, 5] = 40.0
    model.data[0, 6, 5, 5] = 50.0
    model.data[0, 7, 5, 5] = 52.0
    model.data[0, 8, 5, 5] = 63.0
    model.data[0, 9, 5, 5] = 68.0

    for i in range(10, 100):
        model.data[0, i, 5, 5] = i * 5

    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, True)

    assert np.max(out_model.groupdq) == GOOD


def test_saturated_pix(setup_inputs):
    """
    This test is based on an exposure that has some pixels flagged as saturated
    in one or more groups, which the jump step is supposed to ignore, but an
    old version of the code was setting JUMP flags for some of the saturated
    groups. This is to verify that the saturated groups are no longer flagged
    with jumps.
    """
    grouptime = 3.0
    ingain = 1.0
    inreadnoise = 10.7
    ngroups = 7
    nrows = 6
    ncols = 6

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=nrows,
                                       ncols=ncols, gain=ingain,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # Setup the needed input pixel and DQ values
    model.data[0, :, 1, 1] = [639854.75, 4872.451, -17861.791, 14022.15, 22320.176,
                              1116.3828, 1936.9746]
    model.groupdq[0, :, 1, 1] = [0, 0, 0, 0, 0, 2, 2]
    model.data[0, :, 2, 2] = [8.25666812e+05, -1.10471914e+05, 1.95755371e+02, 1.83118457e+03,
                              1.72250879e+03, 1.81733496e+03, 1.65188281e+03]
    model.groupdq[0, :, 2, 2] = [0, 0, 2, 2, 2, 2, 2]
    model.data[0, :, 3, 3] = [1228767., 46392.234, -3245.6553, 7762.413,
                              37190.76, 266611.62, 5072.4434]
    model.groupdq[0, :, 3, 3] = [0, 0, 0, 0, 0, 0, 2]
    model.data[0, :, 4, 4] = [7.5306038e+05, 1.8269953e+04, 1.8352356e+02, 2.1245061e+03,
                              2.0628525e+03, 2.1039399e+03, 2.0069873e+03]
    model.groupdq[0, :, 4, 4] = [0, 0, 2, 2, 2, 2, 2]

    # run jump detection
    out_model = detect_jumps(model, gain, rnoise, rejection_threshold=200.0,
                             max_cores=None, max_jump_to_flag_neighbors=200,
                             min_jump_to_flag_neighbors=10, flag_4_neighbors=True)

    # Check the results. There should not be any pixels with DQ values of 6, which
    # is saturated (2) plus jump (4). All the DQ's should be either just 2 or just 4.
    assert_array_equal(out_model.groupdq[0, :, 1, 1], [0, 4, 4, 4, 0, 2, 2])
    assert_array_equal(out_model.groupdq[0, :, 2, 2], [0, 4, 2, 2, 2, 2, 2])
    assert_array_equal(out_model.groupdq[0, :, 3, 3], [0, 4, 4, 0, 0, 4, 2])
    assert_array_equal(out_model.groupdq[0, :, 4, 4], [0, 4, 2, 2, 2, 2, 2])


@pytest.mark.skip(reason="Test is only used to test performance issue. No need to run every time.")
def test_every_pixel_CR_neighbors_flagged(setup_inputs):
    """"
    A multiprocessing test that has a jump in every pixel. This is used
    to test the performance gain from multiprocessing.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, :, :] = 15.0
    model.data[0, 1, :, :] = 20.0
    model.data[0, 2, :, :] = 25.0
    model.data[0, 3, :, :] = 30.0
    model.data[0, 4, :, :] = 35.0
    model.data[0, 5, :, :] = 140.0
    model.data[0, 6, :, :] = 150.0
    model.data[0, 7, :, :] = 160.0
    model.data[0, 8, :, :] = 170.0
    model.data[0, 9, :, :] = 180.0
    out_model = detect_jumps(model, gain, rnoise, 4.0, 'half', 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 5, 5]) == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 6] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == JUMP_DET


def test_crs_on_edge_with_neighbor_flagging(setup_inputs):
    """"
    A test to make sure that the neighbors of CRs on the edges of the
    array are flagged correctly.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=20, ncols=20,
                                       gain=ingain, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    # CR on 1st row
    model.data[0, 0, 0, 15] = 15.0
    model.data[0, 1, 0, 15] = 20.0
    model.data[0, 2, 0, 15] = 25.0
    model.data[0, 3, 0, 15] = 30.0
    model.data[0, 4, 0, 15] = 35.0
    model.data[0, 5, 0, 15] = 140.0
    model.data[0, 6, 0, 15] = 150.0
    model.data[0, 7, 0, 15] = 160.0
    model.data[0, 8, 0, 15] = 170.0
    model.data[0, 9, 0, 15] = 180.0
    # CR on last row
    model.data[0, 0, -1, 5] = 15.0
    model.data[0, 1, -1, 5] = 20.0
    model.data[0, 2, -1, 5] = 25.0
    model.data[0, 3, -1, 5] = 30.0
    model.data[0, 4, -1, 5] = 35.0
    model.data[0, 5, -1, 5] = 140.0
    model.data[0, 6, -1, 5] = 150.0
    model.data[0, 7, -1, 5] = 160.0
    model.data[0, 8, -1, 5] = 170.0
    model.data[0, 9, -1, 5] = 180.0
    # CR on 1st column
    model.data[0, 0, 5, 0] = 15.0
    model.data[0, 1, 5, 0] = 20.0
    model.data[0, 2, 5, 0] = 25.0
    model.data[0, 3, 5, 0] = 30.0
    model.data[0, 4, 5, 0] = 35.0
    model.data[0, 5, 5, 0] = 140.0
    model.data[0, 6, 5, 0] = 150.0
    model.data[0, 7, 5, 0] = 160.0
    model.data[0, 8, 5, 0] = 170.0
    model.data[0, 9, 5, 0] = 180.0
    # CR on last column
    model.data[0, 0, 15, -1] = 15.0
    model.data[0, 1, 15, -1] = 20.0
    model.data[0, 2, 15, -1] = 25.0
    model.data[0, 3, 15, -1] = 30.0
    model.data[0, 4, 15, -1] = 35.0
    model.data[0, 5, 15, -1] = 140.0
    model.data[0, 6, 15, -1] = 150.0
    model.data[0, 7, 15, -1] = 160.0
    model.data[0, 8, 15, -1] = 170.0
    model.data[0, 9, 15, -1] = 180.0

    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 10, True)

    # flag CR and three neighbors of first row CR
    assert out_model.groupdq[0, 5, 0, 15] == JUMP_DET
    assert out_model.groupdq[0, 5, 1, 15] == JUMP_DET
    assert out_model.groupdq[0, 5, 0, 14] == JUMP_DET
    assert out_model.groupdq[0, 5, 0, 16] == JUMP_DET
    assert out_model.groupdq[0, 5, -1, 15] == 0  # The one not to flag
    # flag CR and three neighbors of last row CR
    assert out_model.groupdq[0, 5, -1, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, -2, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, -1, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, -1, 6] == JUMP_DET
    # flag CR and three neighbors of first column CR
    assert out_model.groupdq[0, 5, 5, 0] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 0] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 0] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 1] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, -1] == 0  # The one not to flag
    # flag CR and three neighbors of last column CR
    assert out_model.groupdq[0, 5, 15, -1] == JUMP_DET
    assert out_model.groupdq[0, 5, 15, -2] == JUMP_DET
    assert out_model.groupdq[0, 5, 16, -1] == JUMP_DET
    assert out_model.groupdq[0, 5, 14, -1] == JUMP_DET


def test_onecr_10_groups(setup_inputs):
    """"
    A test to make sure that neighbors are not flagged when they are not requested to be flagged.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=20, ncols=20,
                                       readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 10, False)

    assert out_model.groupdq[0, 5, 5, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == GOOD
    assert out_model.groupdq[0, 5, 6, 5] == GOOD
    assert out_model.groupdq[0, 5, 5, 6] == GOOD
    assert out_model.groupdq[0, 5, 5, 4] == GOOD


def test_onecr_10_groups_fullarray(setup_inputs):
    """"
    A test that has a cosmic ray in the 5th group for all pixels except column 10. In column
    10 the jump is in the 7th group.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=20, ncols=20,
                                       readnoise=inreadnoise, deltatime=grouptime)

    model.data[0, 0, 5, :] = 15.0
    model.data[0, 1, 5, :] = 20.0
    model.data[0, 2, 5, :] = 25.0
    model.data[0, 3, 5, :] = 30.0
    model.data[0, 4, 5, :] = 35.0
    model.data[0, 5, 5, :] = 140.0
    model.data[0, 6, 5, :] = 150.0
    model.data[0, 7, 5, :] = 160.0
    model.data[0, 8, 5, :] = 170.0
    model.data[0, 9, 5, :] = 180.0
    # move the CR to group 7 for row 10 and make difference be 300
    model.data[0, 3, 5, 10] = 100
    model.data[0, 4, 5, 10] = 130
    model.data[0, 5, 5, 10] = 160
    model.data[0, 6, 5, 10] = 190
    model.data[0, 7, 5, 10] = 400
    model.data[0, 8, 5, 10] = 410
    model.data[0, 9, 5, 10] = 420
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 10, False)

    # The jump is in group 5 for columns 0-9
    assert_array_equal(out_model.groupdq[0, 5, 5, 0:10], JUMP_DET)
    # The jump is in group 7 for column 10
    assert out_model.groupdq[0, 7, 5, 10] == JUMP_DET
    # The jump is in group 5 for columns 11+
    assert_array_equal(out_model.groupdq[0, 5, 5, 11:], JUMP_DET)


def test_onecr_50_groups(setup_inputs):
    """"
    A test with a fifty group integration. There are two jumps in pixel 5,5. One in group 5 and
    one in group 30.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 50
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10, ncols=10,
                                       readnoise=inreadnoise, deltatime=grouptime)

    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    model.data[0, 10:30, 5, 5] = np.arange(190, 290, 5)
    model.data[0, 30:50, 5, 5] = np.arange(500, 600, 5)
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 10, False)

    # CR in group 5
    assert out_model.groupdq[0, 5, 5, 5] == JUMP_DET
    # CR in group 30
    assert out_model.groupdq[0, 30, 5, 5] == JUMP_DET
    # groups in between are not flagged
    assert_array_equal(out_model.groupdq[0, 6:30, 5, 5], GOOD)


def test_single_CR_neighbor_flag(setup_inputs):
    """"
    A single CR in a 10 group exposure. Tests that:
    - if neighbor-flagging is set, the 4 neighboring pixels *ARE* flagged, and
    - if neighbor-flagging is *NOT* set, the 4 neighboring pixels are *NOT* flagged
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=5, ncols=6,
                                       gain=ingain, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 3, 3] = 15.0
    model.data[0, 1, 3, 3] = 20.0
    model.data[0, 2, 3, 3] = 25.0
    model.data[0, 3, 3, 3] = 30.0
    model.data[0, 4, 3, 3] = 35.0
    model.data[0, 5, 3, 3] = 140.0
    model.data[0, 6, 3, 3] = 150.0
    model.data[0, 7, 3, 3] = 160.0
    model.data[0, 8, 3, 3] = 170.0
    model.data[0, 9, 3, 3] = 180.0

    # Flag neighbors
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 3, 3]) == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 3] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 3] == JUMP_DET

    # Do not flag neighbors
    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, False)

    assert np.max(out_model.groupdq[0, 5, 3, 3]) == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 4] == GOOD
    assert out_model.groupdq[0, 5, 3, 2] == GOOD
    assert out_model.groupdq[0, 5, 2, 3] == GOOD
    assert out_model.groupdq[0, 5, 4, 3] == GOOD


def test_adjacent_CRs(setup_inputs):
    """
    Three CRs in a 10 group exposure; the CRs have overlapping neighboring
    pixels. This test makes sure that the correct pixels are flagged.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=15, ncols=6,
                                       gain=ingain, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # Populate arrays for 1st CR, centered at (x=2, y=3)
    x = 2
    y = 3
    model.data[0, 0, y, x] = 15.0
    model.data[0, 1, y, x] = 20.0
    model.data[0, 2, y, x] = 26.0
    model.data[0, 3, y, x] = 30.0
    model.data[0, 4, y, x] = 35.0
    model.data[0, 5, y, x] = 140.0
    model.data[0, 6, y, x] = 150.0
    model.data[0, 7, y, x] = 161.0
    model.data[0, 8, y, x] = 170.0
    model.data[0, 9, y, x] = 180.0

    # Populate arrays for 2nd CR, centered at (x=2, y=2)
    x = 2
    y = 2
    model.data[0, 0, y, x] = 20.0
    model.data[0, 1, y, x] = 30.0
    model.data[0, 2, y, x] = 41.0
    model.data[0, 3, y, x] = 51.0
    model.data[0, 4, y, x] = 62.0
    model.data[0, 5, y, x] = 170.0
    model.data[0, 6, y, x] = 200.0
    model.data[0, 7, y, x] = 231.0
    model.data[0, 8, y, x] = 260.0
    model.data[0, 9, y, x] = 290.0

    # Populate arrays for 3rd CR, centered at (x=3, y=2)
    x = 3
    y = 2
    model.data[0, 0, y, x] = 120.0
    model.data[0, 1, y, x] = 140.0
    model.data[0, 2, y, x] = 161.0
    model.data[0, 3, y, x] = 181.0
    model.data[0, 4, y, x] = 202.0
    model.data[0, 5, y, x] = 70.0
    model.data[0, 6, y, x] = 100.0
    model.data[0, 7, y, x] = 131.0
    model.data[0, 8, y, x] = 160.0
    model.data[0, 9, y, x] = 190.0

    out_model = detect_jumps(model, gain, rnoise, 4.0, 1, 200, 4, True)

    # 1st CR (centered at x=2, y=3)
    assert out_model.groupdq[0, 5, 2, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 1] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 3] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 2] == JUMP_DET

    # 2nd CR (centered at x=2, y=2)
    assert out_model.groupdq[0, 5, 1, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 1] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 3] == JUMP_DET

    # 3rd CR (centered at x=3, y=2)
    assert out_model.groupdq[0, 5, 1, 3] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 4] == JUMP_DET


# Need test for multi-ints near zero with positive and negative slopes

def copy_some_meta_data(objfrom, objto, names):
    # add epydoc copy selected metadata from one datamodel to another
    for n in names:
        if hasattr(objfrom, n):
            v = getattr(objfrom, n)
            setattr(objto, n, v)


@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nrows=102, ncols=103,
               nframes=1, grouptime=1.0, gain=1, deltatime=1):

        # Create a Roman DataModel, recast it as a single integration JWST data model,
        # and modify as needed, return the JWST data model and reference models

        # Create a Roman model with the correct SCI data dimensions
        R_shape = (ngroups, nrows, ncols)
        R_model = R_RampModel(R_shape)

        # Populate the Roman model with some metadata
        R_model.meta.exposure.type = 'WFI_IMAGE'
        R_model.meta.exposure.frame_time = deltatime
        R_model.meta.exposure.ngroups = ngroups
        R_model.meta.exposure.group_time = deltatime
        R_model.meta.exposure.nframes = nframes
        R_model.meta.exposure.groupgap = 0
        R_model.meta.observation.date = Time('2017-11-01T00:00:00')
        R_model.meta.observation.time = Time('2017-11-01T05:13:00')

        # Create a JWST model with the correct SCI data dimensions
        J_shape = (1,) + R_shape # single integration
        J_model = J_RampModel(J_shape)

        # Set mission-specific metadata keywords
        R_model.meta.instrument.name = 'WFI'
        R_model.meta.instrument.detector = 'WFI01'
        J_model.meta.instrument.name = 'NIRCAM'
        J_model.meta.instrument.detector = 'NRCA1'

        # Copy non-mission-specific metadata from Roman model to JWST model,
        #   as the jump code works on JWET models
        mod_from = R_model.meta.exposure
        mod_to = J_model.meta.exposure
        names = ['frame_time','group_time','ngroups','groupgap','nframes']
        copy_some_meta_data(mod_from, mod_to, names)

        mod_from = R_model.meta.observation
        mod_to = J_model.meta.observation
        names = ['date','time']
        copy_some_meta_data(mod_from, mod_to, names)

        # Create gain and readnoise reference models
        gain_model = GainModel((nrows, ncols))
        gain_model.data += gain
        gain_model.meta.instrument.name = 'WFI'

        read_noise_model = ReadNoiseModel((nrows, ncols))
        read_noise_model.data += readnoise
        read_noise_model.meta.instrument.name = 'WFI'

        return J_model, read_noise_model, gain_model

    return _setup
