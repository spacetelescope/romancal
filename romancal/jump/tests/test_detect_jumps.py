import numpy as np
from numpy.testing import assert_array_equal
import pytest

from stcal.jump.jump import detect_jumps

dqflags = {
    "GOOD": 0,
    "DO_NOT_USE": 1,
    "SATURATED": 2,
    "JUMP_DET": 4,
    "NO_GAIN_VALUE": 524288,
}


def test_nocrs_noflux(setup_inputs):
    """
    All pixel values are zero. So slope should be zero
    """
    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=5)

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)

    assert np.max(gdq) == dqflags["GOOD"]


def test_nocrs_noflux_badgain_pixel(setup_inputs):
    """
    All pixel values are zero. So slope should be zero, pixel with bad gain
    should have pixel dq set to 'NO_GAIN_VALUE' and 'DO_NOT_USE'
    """
    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=5)

    gain_2d[7, 7] = -10  # bad gain
    gain_2d[17, 17] = np.nan  # bad gain

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)

    # 2 bits are set for each pixel, so use bitwise_and to check is set
    assert np.bitwise_and(pdq[7, 7], dqflags["NO_GAIN_VALUE"]) == \
        dqflags["NO_GAIN_VALUE"]
    assert np.bitwise_and(pdq[7, 7], dqflags["DO_NOT_USE"]) == \
        dqflags["DO_NOT_USE"]
    assert np.bitwise_and(pdq[17, 17], dqflags["NO_GAIN_VALUE"]) == \
        dqflags["NO_GAIN_VALUE"]
    assert np.bitwise_and(pdq[17, 17], dqflags["DO_NOT_USE"]) == \
        dqflags["DO_NOT_USE"]


def test_onecr_10_groups_neighbors_flagged(setup_inputs):
    """
    A single CR in a 10 group exposure
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10,
                                   ncols=10, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    data[0, 0, 5, 5] = 15.0
    data[0, 1, 5, 5] = 20.0
    data[0, 2, 5, 5] = 25.0
    data[0, 3, 5, 5] = 30.0
    data[0, 4, 5, 5] = 35.0
    data[0, 5, 5, 5] = 140.0
    data[0, 6, 5, 5] = 150.0
    data[0, 7, 5, 5] = 160.0
    data[0, 8, 5, 5] = 170.0
    data[0, 9, 5, 5] = 180.0

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)

    assert np.max(gdq[0, 5, 5, 5]) == dqflags["JUMP_DET"]
    assert gdq[0, 5, 5, 6] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 5, 4] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 6, 5] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 4, 5] == dqflags["JUMP_DET"]


def test_nocr_100_groups_nframes1(setup_inputs):
    """
    NO CR in a 100 group exposure to make sure that frames_per_group is passed
    correctly to twopoint_difference. This test recreates the problem found in
    issue #4571.
    """
    grouptime = 3.0
    ingain = 1
    inreadnoise = 7.0
    ngroups = 100

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10,
                                   ncols=10, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    data[0, 0, 5, 5] = 14.0
    data[0, 1, 5, 5] = 20.0
    data[0, 2, 5, 5] = 27.0
    data[0, 3, 5, 5] = 30.0
    data[0, 4, 5, 5] = 38.0
    data[0, 5, 5, 5] = 40.0
    data[0, 6, 5, 5] = 50.0
    data[0, 7, 5, 5] = 52.0
    data[0, 8, 5, 5] = 63.0
    data[0, 9, 5, 5] = 68.0

    for i in range(10, 100):
        data[0, i, 5, 5] = i * 5

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)

    assert np.max(pdq) == dqflags["GOOD"]


def test_multiple_neighbor_jumps_firstlastbad(setup_inputs):
    """
    This test is based on actual MIRI data that was having the incorrect
    group flagged with JUMP_DET (it was flagging group 2 instead of group 5).
    This makes sure that group 5 is getting flagged.
    Note that the first and last frames/groups are all flagged with
    dqflags["DO_NOT_USE"], due to the application of the first/last frame
    steps.
    """
    grouptime = 3.0
    ingain = 5.5
    inreadnoise = 6.5
    ngroups = 10
    nrows = 10
    ncols = 10

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, err,\
        refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                              ncols=ncols, readnoise=inreadnoise,
                              deltatime=grouptime)

    # Setup the desired pixel values
    data[0, :, 1, 1] = [10019.966, 10057.298, 10078.248, 10096.01, 20241.627,
                        20248.752, 20268.047, 20284.895, 20298.705, 20314.25]
    data[0, :, 1, 2] = [10016.457, 10053.907, 10063.568, 10076.166, 11655.773,
                        11654.063, 11681.795, 11693.763, 11712.788, 11736.994]
    data[0, :, 1, 3] = [10013.259, 10050.348, 10070.398, 10097.658, 10766.534,
                        10787.84,  10802.418, 10818.872, 10832.695, 10861.175]
    data[0, :, 1, 4] = [10016.422, 10053.959, 10070.934, 10090.381, 10104.014,
                        10127.665, 10143.687, 10172.227, 10178.138, 10199.59]
    data[0, :, 2, 1] = [10021.067, 10042.973, 10059.062, 10069.323, 18732.406,
                        18749.602, 18771.908, 18794.695, 18803.223, 18819.523]
    data[0, :, 2, 2] = [10019.651, 10043.371, 10056.423, 10085.121, 40584.703,
                        40606.08, 40619.51, 40629.574, 40641.9, 40660.145]
    data[0, :, 2, 3] = [10021.223, 10042.112, 10052.958, 10067.142, 28188.316,
                        28202.922, 28225.557, 28243.79, 28253.883, 28273.586]
    data[0, :, 2, 4] = [10022.608, 10037.174, 10069.476, 10081.729, 11173.748,
                        11177.344, 11201.127, 11219.607, 11229.468, 11243.174]
    data[0, :, 2, 5] = [10011.095, 10047.422, 10061.066, 10079.375, 10106.405,
                        10116.071, 10129.348, 10136.305, 10161.373, 10181.479]
    data[0, :, 3, 1] = [10011.877, 10052.809, 10075.108, 10085.111, 10397.106,
                        10409.291, 10430.475, 10445.3, 10462.004, 10484.906]
    data[0, :, 3, 2] = [10012.124, 10059.202, 10078.984, 10092.74, 11939.488,
                        11958.45, 11977.5625, 11991.776, 12025.897, 12027.326]
    data[0, :, 3, 3] = [10013.282, 10046.887, 10062.308, 10085.447, 28308.426,
                        28318.957, 28335.55, 28353.832, 28371.746, 28388.848]
    data[0, :, 3, 4] = [10016.784, 10048.249, 10060.097, 10074.606, 21506.082,
                        21522.027, 21542.309, 21558.34, 21576.365, 21595.58]
    data[0, :, 3, 5] = [10014.916, 10052.995, 10063.7705, 10092.866, 10538.075,
                        10558.318, 10570.754, 10597.343, 10608.488, 10628.104]
    data[0, :, 4, 1] = [10017.438, 10038.94, 10057.657, 10069.987, 10090.22,
                        10114.296, 10133.543, 10148.657, 10158.109, 10172.842]
    data[0, :, 4, 2] = [10011.129, 10037.982, 10054.445, 10079.703, 10097.964,
                        10110.593, 10135.701, 10149.448, 10171.771, 10185.874]
    data[0, :, 4, 3] = [10021.109, 10043.658, 10063.909, 10072.364, 10766.232,
                        10774.402, 10790.677, 10809.337, 10833.65, 10849.55]
    data[0, :, 4, 4] = [10023.877, 10035.997, 10052.321, 10077.937, 10529.645,
                        10541.947, 10571.127, 10577.249, 10599.716, 10609.544]

    # Flag first and last frame as DO_NOT_USE
    gdq[0, 0, :, :] = 1
    gdq[0, -1, :, :] = 1

    # run jump detection
    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 200.0, 200.0,
                            200.0, 200, 10, True, dqflags)

    # Verify that the correct groups have been flagged. The entries for pixels
    # 2,2 and 3,3 are the ones that had previously been flagged in group 2
    # instead of group 5.
    assert_array_equal(gdq[0, :, 1, 1], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 1, 2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 1, 3], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 1, 4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 2, 1], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 2, 2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 2, 3], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 2, 4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 3, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 3, 2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 3, 3], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 3, 4], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 4, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 4, 2], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 4, 3], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(gdq[0, :, 4, 4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])


def test_saturated_pix(setup_inputs):
    """
    This test is based on an actual JWST NIRSpec exposure that has some pixels
    flagged as saturated in one or more groups, which the jump step is
    supposed to ignore, but an old version of the code was setting JUMP flags
    for some of the saturated groups. This is to verify that the saturated
    groups are no longer flagged with jumps.
    """
    grouptime = 3.0
    ingain = 1.0
    inreadnoise = 10.7
    ngroups = 7
    nrows = 6
    ncols = 6

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # Setup the needed input pixel and DQ values
    data[0, :, 1, 1] = [639854.75, 4872.451, -17861.791, 14022.15, 22320.176,
                        1116.3828, 1936.9746]
    gdq[0, :, 1, 1] = [0, 0, 0, 0, 0, 2, 2]
    data[0, :, 2, 2] = [8.25666812e+05, -1.10471914e+05, 1.95755371e+02,
                        1.83118457e+03, 1.72250879e+03, 1.81733496e+03,
                        1.65188281e+03]
    gdq[0, :, 2, 2] = [0, 0, 2, 2, 2, 2, 2]
    data[0, :, 3, 3] = [1228767., 46392.234, -3245.6553, 7762.413,
                        37190.76, 266611.62, 5072.4434]
    gdq[0, :, 3, 3] = [0, 0, 0, 0, 0, 0, 2]
    data[0, :, 4, 4] = [7.5306038e+05, 1.8269953e+04, 1.8352356e+02,
                        2.1245061e+03, 2.0628525e+03, 2.1039399e+03,
                        2.0069873e+03]
    gdq[0, :, 4, 4] = [0, 0, 2, 2, 2, 2, 2]

    # run jump detection
    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 200.0, 200.0,
                            200.0, 200, 10, True, dqflags)

    # Check the results. There should not be any pixels with DQ values of 6,
    # which is saturated (2) plus jump (4). All the DQ's should be either
    # just 2 or just 4.
    assert_array_equal(gdq[0, :, 1, 1], [0, 4, 4, 4, 0, 2, 2])
    assert_array_equal(gdq[0, :, 2, 2], [0, 4, 2, 2, 2, 2, 2])
    assert_array_equal(gdq[0, :, 3, 3], [0, 4, 4, 0, 0, 4, 2])
    assert_array_equal(gdq[0, :, 4, 4], [0, 4, 2, 2, 2, 2, 2])


@pytest.mark.skip(reason="Test is only used to test performance issue. No need\
                  to run every time.")
def test_every_pixel_CR_neighbors_flagged(setup_inputs):
    """
    A multiprocessing test that has a jump in every pixel. This is used
    to test the performance gain from multiprocessing.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10,
                                   ncols=10, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    data[0, 0, :, :] = 15.0
    data[0, 1, :, :] = 20.0
    data[0, 2, :, :] = 25.0
    data[0, 3, :, :] = 30.0
    data[0, 4, :, :] = 35.0
    data[0, 5, :, :] = 140.0
    data[0, 6, :, :] = 150.0
    data[0, 7, :, :] = 160.0
    data[0, 8, :, :] = 170.0
    data[0, 9, :, :] = 180.0

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)

    assert np.max(gdq[0, 5, 5, 5]) == dqflags["JUMP_DET"]
    assert gdq[0, 5, 5, 6] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 5, 4] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 6, 5] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 4, 5] == dqflags["JUMP_DET"]


def test_crs_on_edge_with_neighbor_flagging(setup_inputs):
    """
    A test to make sure that the neighbors of CRs on the edges of the
    array are flagged correctly.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    nrows = 20
    ncols = 20

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    # CR on 1st row
    data[0, 0, 0, 15] = 15.0
    data[0, 1, 0, 15] = 20.0
    data[0, 2, 0, 15] = 25.0
    data[0, 3, 0, 15] = 30.0
    data[0, 4, 0, 15] = 35.0
    data[0, 5, 0, 15] = 140.0
    data[0, 6, 0, 15] = 150.0
    data[0, 7, 0, 15] = 160.0
    data[0, 8, 0, 15] = 170.0
    data[0, 9, 0, 15] = 180.0

    # CR on last row
    data[0, 0, -1, 5] = 15.0
    data[0, 1, -1, 5] = 20.0
    data[0, 2, -1, 5] = 25.0
    data[0, 3, -1, 5] = 30.0
    data[0, 4, -1, 5] = 35.0
    data[0, 5, -1, 5] = 140.0
    data[0, 6, -1, 5] = 150.0
    data[0, 7, -1, 5] = 160.0
    data[0, 8, -1, 5] = 170.0
    data[0, 9, -1, 5] = 180.0

    # CR on 1st column
    data[0, 0, 5, 0] = 15.0
    data[0, 1, 5, 0] = 20.0
    data[0, 2, 5, 0] = 25.0
    data[0, 3, 5, 0] = 30.0
    data[0, 4, 5, 0] = 35.0
    data[0, 5, 5, 0] = 140.0
    data[0, 6, 5, 0] = 150.0
    data[0, 7, 5, 0] = 160.0
    data[0, 8, 5, 0] = 170.0
    data[0, 9, 5, 0] = 180.0

    # CR on last column
    data[0, 0, 15, -1] = 15.0
    data[0, 1, 15, -1] = 20.0
    data[0, 2, 15, -1] = 25.0
    data[0, 3, 15, -1] = 30.0
    data[0, 4, 15, -1] = 35.0
    data[0, 5, 15, -1] = 140.0
    data[0, 6, 15, -1] = 150.0
    data[0, 7, 15, -1] = 160.0
    data[0, 8, 15, -1] = 170.0
    data[0, 9, 15, -1] = 180.0

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 10, True, dqflags)

    # flag CR and three neighbors of first row CR
    assert gdq[0, 5, 0, 15] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 1, 15] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 0, 14] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 0, 16] == dqflags["JUMP_DET"]
    assert gdq[0, 5, -1, 15] == 0  # The one not to flag

    # flag CR and three neighbors of last row CR
    assert gdq[0, 5, -1, 5] == dqflags["JUMP_DET"]
    assert gdq[0, 5, -2, 5] == dqflags["JUMP_DET"]
    assert gdq[0, 5, -1, 4] == dqflags["JUMP_DET"]
    assert gdq[0, 5, -1, 6] == dqflags["JUMP_DET"]

    # flag CR and three neighbors of first column CR
    assert gdq[0, 5, 5, 0] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 6, 0] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 4, 0] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 5, 1] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 5, -1] == 0  # The one not to flag

    # flag CR and three neighbors of last column CR
    assert gdq[0, 5, 15, -1] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 15, -2] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 16, -1] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 14, -1] == dqflags["JUMP_DET"]


def test_onecr_10_groups(setup_inputs):
    """
    A test to make sure that neighbors are not flagged when they are not
    requested to be flagged.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    nrows = 20
    ncols = 20

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    data[0, 0, 5, 5] = 15.0
    data[0, 1, 5, 5] = 20.0
    data[0, 2, 5, 5] = 25.0
    data[0, 3, 5, 5] = 30.0
    data[0, 4, 5, 5] = 35.0
    data[0, 5, 5, 5] = 140.0
    data[0, 6, 5, 5] = 150.0
    data[0, 7, 5, 5] = 160.0
    data[0, 8, 5, 5] = 170.0
    data[0, 9, 5, 5] = 180.0

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 10, False, dqflags)

    assert gdq[0, 5, 5, 5] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 4, 5] == dqflags["GOOD"]
    assert gdq[0, 5, 6, 5] == dqflags["GOOD"]
    assert gdq[0, 5, 5, 6] == dqflags["GOOD"]
    assert gdq[0, 5, 5, 4] == dqflags["GOOD"]


def test_onecr_10_groups_fullarray(setup_inputs):
    """
    A test that has a cosmic ray in the 5th group for all pixels except column
    10. In column 10 the jump is in the 7th group.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    nrows = 20
    ncols = 20

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    data[0, 0, 5, :] = 15.0
    data[0, 1, 5, :] = 20.0
    data[0, 2, 5, :] = 25.0
    data[0, 3, 5, :] = 30.0
    data[0, 4, 5, :] = 35.0
    data[0, 5, 5, :] = 140.0
    data[0, 6, 5, :] = 150.0
    data[0, 7, 5, :] = 160.0
    data[0, 8, 5, :] = 170.0
    data[0, 9, 5, :] = 180.0

    # move the CR to group 7 for row 10 and make difference be 300
    data[0, 3, 5, 10] = 100
    data[0, 4, 5, 10] = 130
    data[0, 5, 5, 10] = 160
    data[0, 6, 5, 10] = 190
    data[0, 7, 5, 10] = 400
    data[0, 8, 5, 10] = 410
    data[0, 9, 5, 10] = 420

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 10, False, dqflags)

    # The jump is in group 5 for columns 0-9
    assert_array_equal(gdq[0, 5, 5, 0:10], dqflags["JUMP_DET"])

    # The jump is in group 7 for column 10
    assert gdq[0, 7, 5, 10] == dqflags["JUMP_DET"]

    # The jump is in group 5 for columns 11+
    assert_array_equal(gdq[0, 5, 5, 11:], dqflags["JUMP_DET"])


def test_onecr_50_groups(setup_inputs):
    """
    A test with a fifty group integration. There are two jumps in pixel 5,5.
    One in group 5 and one in group 30.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 50
    nrows = 10
    ncols = 10

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    data[0, 0, 5, 5] = 15.0
    data[0, 1, 5, 5] = 20.0
    data[0, 2, 5, 5] = 25.0
    data[0, 3, 5, 5] = 30.0
    data[0, 4, 5, 5] = 35.0
    data[0, 5, 5, 5] = 140.0
    data[0, 6, 5, 5] = 150.0
    data[0, 7, 5, 5] = 160.0
    data[0, 8, 5, 5] = 170.0
    data[0, 9, 5, 5] = 180.0
    data[0, 10:30, 5, 5] = np.arange(190, 290, 5)
    data[0, 30:50, 5, 5] = np.arange(500, 600, 5)

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 10, False, dqflags)

    # CR in group 5
    assert gdq[0, 5, 5, 5] == dqflags["JUMP_DET"]

    # CR in group 30
    assert gdq[0, 30, 5, 5] == dqflags["JUMP_DET"]

    # groups in between are not flagged
    assert_array_equal(gdq[0, 6:30, 5, 5], dqflags["GOOD"])


def test_single_CR_neighbor_flag(setup_inputs):
    """
    A single CR in a 10 group exposure. Tests that:
    - if neighbor-flagging is set, the 4 neighboring pixels *ARE* flagged, and
    - if neighbor-flagging is *NOT* set, the 4 neighboring pixels are *NOT*
      flagged
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    nrows = 5
    ncols = 6

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    data[0, 0, 3, 3] = 15.0
    data[0, 1, 3, 3] = 20.0
    data[0, 2, 3, 3] = 25.0
    data[0, 3, 3, 3] = 30.0
    data[0, 4, 3, 3] = 35.0
    data[0, 5, 3, 3] = 140.0
    data[0, 6, 3, 3] = 150.0
    data[0, 7, 3, 3] = 160.0
    data[0, 8, 3, 3] = 170.0
    data[0, 9, 3, 3] = 180.0

    # Flag neighbors
    gdq_1, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                              gain_2d, rn_2d, 4.0, 5.0,
                              6.0, 200, 4, True, dqflags)

    assert np.max(gdq_1[0, 5, 3, 3]) == dqflags["JUMP_DET"]
    assert gdq_1[0, 5, 3, 4] == dqflags["JUMP_DET"]
    assert gdq_1[0, 5, 3, 2] == dqflags["JUMP_DET"]
    assert gdq_1[0, 5, 2, 3] == dqflags["JUMP_DET"]
    assert gdq_1[0, 5, 4, 3] == dqflags["JUMP_DET"]

    # Do not flag neighbors
    gdq_2, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                              gain_2d, rn_2d, 4.0, 5.0,
                              6.0, 200, 4, False, dqflags)

    assert np.max(gdq_2[0, 5, 3, 3]) == dqflags["JUMP_DET"]
    assert gdq_2[0, 5, 3, 4] == dqflags["GOOD"]
    assert gdq_2[0, 5, 3, 2] == dqflags["GOOD"]
    assert gdq_2[0, 5, 2, 3] == dqflags["GOOD"]
    assert gdq_2[0, 5, 4, 3] == dqflags["GOOD"]


'''
#   skip next as it's for multiprocess
def SKIP_test_proc(setup_inputs):
    """
    A single CR in a 10 group exposure. Verify that the pixels flagged using
    multiprocessing are identical to the pixels flagged when no
    multiprocessing is done.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10,
                                   ncols=10, readnoise=inreadnoise,
                                   deltatime=grouptime)

    data[0, 0, 2, 3] = 15.0
    data[0, 1, 2, 3] = 21.0
    data[0, 2, 2, 3] = 25.0
    data[0, 3, 2, 3] = 30.2
    data[0, 4, 2, 3] = 35.0
    data[0, 5, 2, 3] = 140.0
    data[0, 6, 2, 3] = 151.0
    data[0, 7, 2, 3] = 160.0
    data[0, 8, 2, 3] = 170.0
    data[0, 9, 2, 3] = 180.0

    gdq_a, pdq_a = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)



    out_model_b = detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 200, 4,
                               True, dqflags)
    assert_array_equal(out_model_a.groupdq, out_model_b.groupdq)

    out_model_c = detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 200, 4,
                               True, dqflags)
    assert_array_equal(out_model_a.groupdq, out_model_c.groupdq)
'''

def test_adjacent_CRs(setup_inputs):
    """
    Three CRs in a 10 group exposure; the CRs have overlapping neighboring
    pixels. This test makes sure that the correct pixels are flagged.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    nrows = 15
    ncols = 6

    frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, pdq, \
        err, refout = setup_inputs(ngroups=ngroups, gain=ingain, nrows=nrows,
                                   ncols=ncols, readnoise=inreadnoise,
                                   deltatime=grouptime)

    # Populate arrays for 1st CR, centered at (x=2, y=3)
    x = 2
    y = 3
    data[0, 0, y, x] = 15.0
    data[0, 1, y, x] = 20.0
    data[0, 2, y, x] = 26.0
    data[0, 3, y, x] = 30.0
    data[0, 4, y, x] = 35.0
    data[0, 5, y, x] = 140.0
    data[0, 6, y, x] = 150.0
    data[0, 7, y, x] = 161.0
    data[0, 8, y, x] = 170.0
    data[0, 9, y, x] = 180.0

    # Populate arrays for 2nd CR, centered at (x=2, y=2)
    x = 2
    y = 2
    data[0, 0, y, x] = 20.0
    data[0, 1, y, x] = 30.0
    data[0, 2, y, x] = 41.0
    data[0, 3, y, x] = 51.0
    data[0, 4, y, x] = 62.0
    data[0, 5, y, x] = 170.0
    data[0, 6, y, x] = 200.0
    data[0, 7, y, x] = 231.0
    data[0, 8, y, x] = 260.0
    data[0, 9, y, x] = 290.0

    # Populate arrays for 3rd CR, centered at (x=3, y=2)
    x = 3
    y = 2
    data[0, 0, y, x] = 120.0
    data[0, 1, y, x] = 140.0
    data[0, 2, y, x] = 161.0
    data[0, 3, y, x] = 181.0
    data[0, 4, y, x] = 202.0
    data[0, 5, y, x] = 70.0
    data[0, 6, y, x] = 100.0
    data[0, 7, y, x] = 131.0
    data[0, 8, y, x] = 160.0
    data[0, 9, y, x] = 190.0

    gdq, pdq = detect_jumps(frames_per_group, data, gdq, pdq, err,
                            gain_2d, rn_2d, 4.0, 5.0,
                            6.0, 200, 4, True, dqflags)

    # 1st CR (centered at x=2, y=3)
    assert gdq[0, 5, 2, 2] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 3, 1] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 3, 2] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 3, 3] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 4, 2] == dqflags["JUMP_DET"]

    # 2nd CR (centered at x=2, y=2)
    assert gdq[0, 5, 1, 2] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 2, 1] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 2, 3] == dqflags["JUMP_DET"]

    # 3rd CR (centered at x=3, y=2)
    assert gdq[0, 5, 1, 3] == dqflags["JUMP_DET"]
    assert gdq[0, 5, 2, 4] == dqflags["JUMP_DET"]


# Need test for multi-ints near zero with positive and negative slopes

@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nints=1, nrows=102, ncols=103,
               frames_per_group=1, grouptime=1.0, gain=1, deltatime=1,
               subarray=False):

        if subarray:
            shape = (nints, ngroups, 20, 20)
        else:
            shape = (nints, ngroups, nrows, ncols)

        refout_shape = (nints, ngroups, nrows, int(ncols/4))

        ncols = shape[3]  # = ncols, xsize
        nrows = shape[2]   # = nrows, ysize
        gain_2d = np.zeros((nrows, ncols)) + gain
        rn_2d = np.zeros((nrows, ncols)) + readnoise

        data = np.zeros(shape, dtype=np.float)
        gdq = np.zeros(shape, dtype=np.int32)
        pdq = np.zeros((nrows, ncols), dtype=np.int32)

        err = np.zeros(shape, dtype=np.float)

        refout = np.zeros(refout_shape, dtype=np.float)

        return frames_per_group, ncols, nrows, gain_2d, rn_2d, data, gdq, \
            pdq, err, refout

    return _setup
