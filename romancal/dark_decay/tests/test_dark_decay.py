import numpy as np
from roman_datamodels import datamodels

from romancal.dark_decay.dark_decay import subtract_dark_decay
from romancal.dark_decay.dark_decay_step import DarkDecayStep


def create_ramp_model(nresultants, nrows=4096, ncols=4096):
    # make a ramp model with fake data
    ramp_model = datamodels.RampModel.create_fake_data(
        shape=(nresultants, nrows, ncols)
    )

    # add necessary metadata
    ramp_model.meta.exposure.frame_time = 1.0

    # add required cal steps
    ramp_model.meta.cal_step = {}
    for step_name in ramp_model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        ramp_model.meta.cal_step[step_name] = "INCOMPLETE"

    # add a read pattern with an increasing number of groups
    read_pattern = []
    next_value = 1
    for n in range(nresultants):
        read_pattern.append(list(range(next_value, next_value + n + 1)))
        next_value = read_pattern[-1][-1] + 1

    ramp_model.meta.exposure.read_pattern = read_pattern
    ramp_model.pixeldq = np.zeros((nrows, ncols), dtype=ramp_model.pixeldq.dtype)
    ramp_model.groupdq = np.zeros(
        (nresultants, nrows, ncols), dtype=ramp_model.groupdq.dtype
    )

    return ramp_model


def test_dark_decay():
    model = create_ramp_model(2, nrows=4096, ncols=4096)
    model.meta.instrument.detector = "WFI01"
    origdata = model.data.copy()
    decayref = datamodels.DarkdecaysignalRefModel.create_fake_data()
    amplitude = 1
    time_constant = 23
    decayref.decay_table["WFI01"] = dict(
        amplitude=amplitude, time_constant=time_constant
    )
    DarkDecayStep.call(model, override_darkdecaysignal=decayref)
    delta = model.data - origdata

    # check that we did just subtract off the relevant signal;
    # i.e., the calling the step does what we expect
    expectation = np.zeros_like(model.data)
    subtract_dark_decay(
        expectation,
        amplitude,
        time_constant,
        model.meta.exposure.frame_time,
        model.meta.exposure.read_pattern,
        1,
    )
    assert np.allclose(expectation, delta)

    # check that if the amplitude doubles the subtraction doubles
    expectation = np.zeros_like(model.data)
    subtract_dark_decay(
        expectation,
        2 * amplitude,
        time_constant,
        model.meta.exposure.frame_time,
        model.meta.exposure.read_pattern,
        1,
    )
    assert np.allclose(expectation, delta * 2)

    # check that early times have larger corrections than late times
    # (i.e., the signal decays)
    assert np.all(np.abs(expectation[0]) >= np.abs(expectation[1]))

    # check that the signal dilutes when adding later frames to a resultant
    # expectation1: resultant 0 is frame 1 only
    # expectation2: resultant 0 is average of frames 1 and 2
    # averaging with a later frame should reduce the dark decay signal
    expectation1 = np.zeros_like(model.data)
    subtract_dark_decay(
        expectation1,
        amplitude,
        time_constant,
        model.meta.exposure.frame_time,
        [[1], [2, 3]],
        1,
    )
    expectation2 = np.zeros_like(model.data)
    subtract_dark_decay(
        expectation2,
        amplitude,
        time_constant,
        model.meta.exposure.frame_time,
        [[1, 2], [3]],
        1,
    )
    assert np.all(np.abs(expectation1[0]) >= np.abs(expectation2[0]))

    # check that the amplitude is about right in the first read.
    # frame time of 1 is much smaller than the time constant of 23.
    # amplitude should correspond closely to an average of the first read
    assert np.abs(np.mean(expectation1[0]) + amplitude) < 0.001
