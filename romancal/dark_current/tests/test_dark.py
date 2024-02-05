"""
Unit tests for dark current correction
"""

import numpy as np
import pytest
import roman_datamodels as rdm
from astropy import units as u
from roman_datamodels import maker_utils
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
    assert type(result) == RampModel
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
        ramp_model.data[0, 0, i] = i * ramp_model.data.unit
        darkref_model.data[0, 0, i] = i * 0.1 * darkref_model.data.unit
    orig_model = ramp_model.copy()

    # Perform Dark Current subtraction step
    result = DarkCurrentStep.call(ramp_model, override_dark=darkref_model)

    # check that the dark file is subtracted frame by frame from the science data
    diff = orig_model.data.value - darkref_model.data.value

    # test that the output data file is equal to the difference found when subtracting
    # reffile from sci file
    np.testing.assert_array_equal(
        result.data.value, diff, err_msg="dark file should be subtracted from sci file "
    )


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_output_dark_file(tmpdir, instrument, exptype):
    """Test that the the step can output a proper (optional) dark file"""
    path = str(tmpdir / "dark_out.asdf")

    # Set test size
    shape = (2, 20, 20)

    # Create test ramp and dark models
    ramp_model, darkref_model = create_ramp_and_dark(shape, instrument, exptype)

    # Perform Dark Current subtraction step
    DarkCurrentStep.call(ramp_model, override_dark=darkref_model, dark_output=path)

    # Open dark file
    dark_out_file_model = rdm.open(path)

    # Test dark file results
    assert type(dark_out_file_model) == DarkRefModel
    assert dark_out_file_model.validate() is None
    assert dark_out_file_model.data.shape == shape
    assert dark_out_file_model.dq.shape == shape[1:]


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_getbestrefs(tmpdir, instrument, exptype):
    """Test that the the step will skip if CRDS returns N/A for the ref file"""
    path = str(tmpdir / "dark_out.asdf")

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
    assert type(dark_out_file_model) == DarkRefModel
    assert dark_out_file_model.validate() is None
    assert dark_out_file_model.data.shape == shape
    assert dark_out_file_model.dq.shape == shape[1:]


def create_ramp_and_dark(shape, instrument, exptype):
    """Helper function to create test ramp and dark models"""

    # Create test ramp model
    ramp = maker_utils.mk_ramp(shape=shape)
    ramp.meta.instrument.name = instrument
    ramp.meta.instrument.detector = "WFI01"
    ramp.meta.instrument.optical_element = "F158"
    ramp.meta.exposure.type = exptype
    ramp.data = u.Quantity(np.ones(shape, dtype=np.float32), u.DN, dtype=np.float32)
    ramp_model = RampModel(ramp)

    # Create dark model
    darkref = maker_utils.mk_dark(shape=shape)
    darkref_model = DarkRefModel(darkref)

    return ramp_model, darkref_model
