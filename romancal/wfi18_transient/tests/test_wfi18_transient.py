import numpy as np
from roman_datamodels import datamodels, dqflags

from romancal.wfi18_transient.wfi18_transient import _double_exp, _frame_read_times
from romancal.wfi18_transient.wfi18_transient_step import WFI18TransientStep

MASK_FLAG = dqflags.group.DO_NOT_USE | dqflags.group.WFI18_TRANSIENT


def create_ramp_model(nresultants, nrows=4096, ncols=4096):
    # make a ramp model with fake data
    ramp_model = datamodels.RampModel.create_fake_data(
        shape=(nresultants, nrows, ncols)
    )

    # add necessary metadata
    ramp_model.meta.instrument.detector = "WFI18"
    ramp_model.meta.exposure.frame_time = 1.0
    ramp_model.meta.exposure.sca_number = 18

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


def transient_glow():
    read_times = _frame_read_times(18, 1.0)
    glow = _double_exp(read_times, 0.1, 0.01, 0.01, 0.01)
    return glow


def test_wfi18_transient(caplog):
    model = create_ramp_model(5)

    # Fill the data array with a flat value
    model.data[:] = 1.0

    # Add a glow to the bottom of the detector in the first read,
    # not including the reference pixels
    model.data[0, 4:-4, 4:-4] += transient_glow()[4:-4, 4:-4]
    assert not np.allclose(model.data, 1.0, atol=1e-5)

    # Correct out the glow
    result = WFI18TransientStep.call(model)
    assert result.meta.cal_step.wfi18_transient == "COMPLETE"

    # For this test data, the fit should succeed. The correction
    # should be pretty good, restoring the value to the expected
    # value of 1.0
    np.testing.assert_allclose(result.data, 1.0, atol=1e-5)
    np.testing.assert_equal(result.groupdq, 0)


def test_wfi18_transient_flat_data(caplog):
    model = create_ramp_model(5)

    # Fill the data array with a flat value
    model.data[:] = 1.0

    result = WFI18TransientStep.call(model)
    assert result.meta.cal_step.wfi18_transient == "COMPLETE"

    # For flat data, the fit should succeed and there
    # should be no change in the output
    np.testing.assert_equal(result.data, 1)
    np.testing.assert_equal(result.groupdq, 0)


def test_wfi18_transient_wrong_detector():
    model = create_ramp_model(1)
    model.meta.instrument.detector = "WFI01"

    # attempting to call on the wrong detector will skip the step
    # and set the cal step status to N/A
    result = WFI18TransientStep.call(model)
    assert result.meta.cal_step.wfi18_transient == "N/A"


def test_wfi18_transient_too_few_resultants(caplog):
    model = create_ramp_model(1)

    # attempting to call with resultants < 5 will fall back on masking
    result = WFI18TransientStep.call(model)
    assert result.meta.cal_step.wfi18_transient == "COMPLETE"

    assert "Masking affected rows" in caplog.text
    np.testing.assert_equal(result.data, 0)
    np.testing.assert_equal(result.groupdq[0, :1000], MASK_FLAG)
    np.testing.assert_equal(result.groupdq[0, 1000:], 0)
    np.testing.assert_equal(result.groupdq[1:], 0)


def test_wfi18_transient_fit_failure(caplog):
    model = create_ramp_model(5)
    model.data[:] = np.nan

    # attempting to call with invalid data will fall back on masking
    result = WFI18TransientStep.call(model)
    assert result.meta.cal_step.wfi18_transient == "COMPLETE"

    assert "Transient fit failed; masking affected rows instead" in caplog.text
    np.testing.assert_equal(result.data, np.nan)
    np.testing.assert_equal(result.groupdq[0, :1000], MASK_FLAG)
    np.testing.assert_equal(result.groupdq[0, 1000:], 0)
    np.testing.assert_equal(result.groupdq[1:], 0)


def test_wfi18_transient_save_results(tmp_path):
    # Make a small model
    model = create_ramp_model(1, nrows=10, ncols=10)
    model.meta.filename = "test_refpix.asdf"
    input_filename = str(tmp_path / model.meta.filename)
    model.save(input_filename)

    # Run on the file on disk, save the output
    WFI18TransientStep.call(input_filename, save_results=True, output_dir=str(tmp_path))

    # Output should have the correct suffix
    expected_output = tmp_path / "test_wfi18_transient.asdf"
    assert expected_output.exists()

    # Saved output has the expected type and metadata
    with datamodels.open(expected_output) as output_model:
        assert isinstance(output_model, datamodels.RampModel)
        assert output_model.meta.cal_step.wfi18_transient == "COMPLETE"
