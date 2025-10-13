"""Unit tests for dark current correction"""

import numpy as np
import pytest
import roman_datamodels as rdm
from roman_datamodels.datamodels import DarkRefModel, ImageModel

from romancal.dark_current import DarkCurrentStep


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a DARK reffile"""
    # Set test size, data + reference files
    ref_shape = (10, 28, 28)
    image_shape = (ref_shape[1], ref_shape[2])

    # Create test rampfit and dark models
    rampfit_model, darkref_model = create_image_and_dark(ref_shape, instrument, exptype)
    # Perform Dark Current subtraction step
    result = DarkCurrentStep.call(rampfit_model, override_dark=darkref_model)

    # Test dark results
    trim_shape = (image_shape[0] - 8, image_shape[1] - 8)  # size of the trimmed image
    assert (result.data == rampfit_model.data).all()
    assert isinstance(result, ImageModel)
    assert result.validate() is None
    assert result.data.shape == trim_shape
    assert result.dq.shape == trim_shape
    assert result.meta.cal_step.dark == "COMPLETE"
    assert result.data.dtype == np.float32
    assert result.dq.dtype == np.uint32


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dark_step_subtraction(instrument, exptype):
    """Test that the values in a dark reference file are properly subtracted"""
    # Set test size
    ref_shape = (10, 28, 28)

    # Create test ramp and dark models
    image_model, darkref_model = create_image_and_dark(ref_shape, instrument, exptype)

    # populate data array of science cube
    orig_model = image_model.copy()

    # Perform Dark Current subtraction step
    result = DarkCurrentStep.call(image_model, override_dark=darkref_model)

    diff = orig_model.data - (darkref_model.dark_slope[4:-4, 4:-4])

    # test that the output data file is equal to the difference found when subtracting
    # reffile from sci file
    assert np.all(np.isclose(result.data, diff, atol=1.0e-7))


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
    ref_shape = (10, 28, 28)
    image_shape = (ref_shape[1], ref_shape[2])

    # Create test ramp and dark models
    image_model, darkref_model = create_image_and_dark(ref_shape, instrument, exptype)

    # Perform Dark Current subtraction step
    DarkCurrentStep.call(image_model, override_dark=darkref_model, dark_output=path)

    # Open dark file
    dark_out_file_model = rdm.open(path)

    # Test dark file results
    assert isinstance(dark_out_file_model, DarkRefModel)
    assert dark_out_file_model.validate() is None
    assert dark_out_file_model.dq.shape == image_shape


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
    ref_shape = (10, 28, 28)
    image_shape = (ref_shape[1], ref_shape[2])

    # Create test ramp and dark models
    image_model, darkref_model = create_image_and_dark(ref_shape, instrument, exptype)

    # Perform Dark Current subtraction step with override = N/A
    result = DarkCurrentStep.call(image_model, override_dark="N/A")
    assert result.meta.cal_step.dark == "SKIPPED"

    # Perform Dark Current subtraction step
    DarkCurrentStep.call(image_model, override_dark=darkref_model, dark_output=path)

    # Open dark file
    dark_out_file_model = rdm.open(path)

    # Test dark file results
    assert isinstance(dark_out_file_model, DarkRefModel)
    assert dark_out_file_model.validate() is None
    assert dark_out_file_model.dq.shape == image_shape


def create_image_and_dark(shape, instrument, exptype):
    """Helper function to create test image and dark models"""
    # Create test image model
    image_shape = (shape[1], shape[2])
    image = ImageModel.create_fake_data(shape=image_shape)
    image.meta.cal_step = {}
    for step_name in image.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        image.meta.cal_step[step_name] = "INCOMPLETE"
    image.meta.cal_logs = []
    image.meta.instrument.name = instrument
    image.meta.instrument.detector = "WFI01"
    image.meta.instrument.optical_element = "F158"
    image.meta.exposure.type = exptype
    image.meta.exposure.read_pattern = [[1], [2, 3], [4], [5, 6, 7, 8], [9, 10], [11]]
    image.data = np.ones(image_shape, dtype=np.float32)
    image.data = np.full(image_shape, 0.2, dtype=np.float32)
    image.dq = np.zeros(image_shape, dtype=image.dq.dtype)
    # Trim the image
    image.data = image.data[4:-4, 4:-4]
    image.dq = image.dq[4:-4, 4:-4]
    image.err = image.err[4:-4, 4:-4]

    # Create dark model
    darkref = DarkRefModel.create_fake_data(shape=shape)
    darkref.dark_slope = np.full((image_shape), 5.3e-03, dtype=np.float32)
    darkref.dark_slope_error = np.full((image_shape), 2.6e-05, dtype=np.float32)
    darkref.dq = np.zeros(image_shape, dtype=image.dq.dtype)

    return image, darkref
