"""
Unit tests for dark current correction
"""

import numpy as np
import pytest
import roman_datamodels as rdm
from roman_datamodels import stnode as st
from roman_datamodels.datamodels import DarkRefModel, RampModel

from romancal.dark_current import DarkCurrentStep


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a DARK reffile"""

    # Set test size
    shape = (2, 20, 20)

    # Create test ramp and dark models
    ramp_model, darkref_model = create_ramp_and_dark(shape, instrument, exptype)

    # Perform Dark Current subtraction step
    result = DarkCurrentStep.call(ramp_model, override_dark=darkref_model)

    # Test dark results
    assert (result.data == ramp_model.data).all()
    assert isinstance(result, RampModel)
    assert result.validate() is None
    assert result.data.shape == shape
    assert result.groupdq.shape == shape
    assert result.pixeldq.shape == shape[1:]
    assert result.meta.cal_step.dark == "COMPLETE"
    assert result.data.dtype == np.float32
    assert result.pixeldq.dtype == np.uint32
    assert result.groupdq.dtype == np.uint8


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_subtraction(instrument, exptype):
    """Test that the values in a dark reference file are properly subtracted"""

    # Set test size
    shape = (2, 20, 20)

    # Create test ramp and dark models
    ramp_model, darkref_model = create_ramp_and_dark(shape, instrument, exptype)

    # populate data array of science cube
    for i in range(0, 20):
        ramp_model.data[0, 0, i] = i
        darkref_model.data[0, 0, i] = i * 0.1
    orig_model = ramp_model.copy()

    # Perform Dark Current subtraction step
    result = DarkCurrentStep.call(ramp_model, override_dark=darkref_model)

    # check that the dark file is subtracted frame by frame from the science data
    diff = orig_model.data - darkref_model.data

    # test that the output data file is equal to the difference found when subtracting
    # reffile from sci file
    np.testing.assert_array_equal(
        result.data, diff, err_msg="dark file should be subtracted from sci file "
    )


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_output_dark_file(tmp_path, instrument, exptype):
    """Test that the the step can output a proper (optional) dark file"""
    path = str(tmp_path / "dark_out.asdf")

    # Set test size
    shape = (2, 20, 20)

    # Create test ramp and dark models
    ramp_model, darkref_model = create_ramp_and_dark(shape, instrument, exptype)

    # Perform Dark Current subtraction step
    DarkCurrentStep.call(ramp_model, override_dark=darkref_model, dark_output=path)

    # Open dark file
    dark_out_file_model = rdm.open(path)

    # Test dark file results
    assert isinstance(dark_out_file_model, DarkRefModel)
    assert dark_out_file_model.validate() is None
    assert dark_out_file_model.data.shape == shape
    assert dark_out_file_model.dq.shape == shape[1:]


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_getbestrefs(tmp_path, instrument, exptype):
    """Test that the the step will skip if CRDS returns N/A for the ref file"""
    path = str(tmp_path / "dark_out.asdf")

    # Set test size
    shape = (2, 20, 20)

    # Create test ramp and dark models
    ramp_model, darkref_model = create_ramp_and_dark(shape, instrument, exptype)

    # Perform Dark Current subtraction step with override = N/A
    result = DarkCurrentStep.call(ramp_model, override_dark="N/A")
    assert result.meta.cal_step.dark == "SKIPPED"

    # Perform Dark Current subtraction step
    DarkCurrentStep.call(ramp_model, override_dark=darkref_model, dark_output=path)

    # Open dark file
    dark_out_file_model = rdm.open(path)

    # Test dark file results
    assert isinstance(dark_out_file_model, DarkRefModel)
    assert dark_out_file_model.validate() is None
    assert dark_out_file_model.data.shape == shape
    assert dark_out_file_model.dq.shape == shape[1:]


def create_ramp_and_dark(shape, instrument, exptype):
    """Helper function to create test ramp and dark models"""

    # Create test ramp model
    ramp = RampModel.create_fake_data(shape=shape)
    ramp.meta.cal_step = st.L2CalStep.create_fake_data()
    ramp.meta.instrument.name = instrument
    ramp.meta.instrument.detector = "WFI01"
    ramp.meta.instrument.optical_element = "F158"
    ramp.meta.exposure.type = exptype
    ramp.meta.exposure.read_pattern = [[1], [2, 3], [4], [5, 6, 7, 8], [9, 10], [11]]
    ramp.data = np.ones(shape, dtype=np.float32)
    ramp.pixeldq = np.zeros(shape[1:], dtype=ramp.pixeldq.dtype)

    # Create dark model
    darkref = DarkRefModel.create_fake_data(shape=shape[1:])
    darkref.data = np.zeros(shape, dtype=darkref.data.dtype)

    return ramp, darkref
