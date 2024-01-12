""" Roman tests for flat field correction """

import pytest
import roman_datamodels as rdm
from metrics_logger.decorators import metrics_logger

from romancal.pipeline.exposure_pipeline import ExposurePipeline


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS413")
def test_16resultants_image_processing(rtdata, ignore_asdf_paths):
    """Tests for imaging processing requirements for 16 resultants (DMS413)"""
    # The input data is from INS for stress testing at some point this should be generated
    # every time new data is needed.

    input_dark = "roman_dark_WFI01_IMAGE_STRESS_TEST_16_MA_TABLE_998_D1.asdf"
    rtdata.get_data(f"WFI/image/{input_dark}")

    input_data = "r00r1601001001001001_01101_0001_WFI01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r00r1601001001001001_01101_0001_WFI01_cal.asdf"
    rtdata.output = output
    args = [
        "--disable-crds-steppars",
        "--steps.dark_current.override_dark=roman_dark_WFI01_IMAGE_STRESS_TEST_16_MA_TABLE_998_D1.asdf",
        "roman_elp",
        rtdata.input,
    ]
    ExposurePipeline.from_cmdline(args)

    # Perform DMS tests
    # Initial prep
    pipeline = ExposurePipeline()
    # rtdata.get_data(f"WFI/image/{output}")
    model = rdm.open(output, lazy_load=False)
    uncal_data = rdm.open(input_data, lazy_load=False)

    # DMS280 result is an ImageModel
    pipeline.log.info(
        "DMS413 MSG: Testing that result is a Level 2 model......."
        + passfail(isinstance(model, rdm.datamodels.ImageModel))
    )

    pipeline.log.info(
        "DMS413 MSG: Testing that there are 16 resultants in the input file......."
        + passfail(len(uncal_data.meta.exposure.read_pattern) == 16)
    )

    pipeline.log.info(
        "DMS413 MSG: Testing that there are 16 resultants listed in the output file......."
        + passfail(len(model.meta.exposure.read_pattern) == 16)
    )

    # Ensure step completion is as expected
    assert model.meta.cal_step.dq_init == "COMPLETE"
    assert model.meta.cal_step.saturation == "COMPLETE"
    assert model.meta.cal_step.linearity == "COMPLETE"
    assert model.meta.cal_step.dark == "COMPLETE"
    assert model.meta.cal_step.jump == "COMPLETE"
    assert model.meta.cal_step.ramp_fit == "COMPLETE"
    assert model.meta.cal_step.assign_wcs == "COMPLETE"
    assert model.meta.cal_step.flat_field == "COMPLETE"
    assert model.meta.cal_step.photom == "COMPLETE"


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS414")
def test_16resultants_spectral_processing(rtdata, ignore_asdf_paths):
    """Tests for imaging processing requirements for 16 resultants (DMS413)"""
    # The input data is from INS for stress testing at some point this should be generated
    # by INS every time new data is needed.

    input_dark = "roman_dark_WFI01_IMAGE_STRESS_TEST_16_MA_TABLE_998_D1.asdf"
    rtdata.get_data(f"WFI/image/{input_dark}")

    input_data = "r10r1601001001001001_01101_0001_WFI01_uncal.asdf"
    rtdata.get_data(f"WFI/grism/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r10r1601001001001001_01101_0001_WFI01_cal.asdf"
    rtdata.output = output
    args = [
        "--disable-crds-steppars",
        "--steps.dark_current.override_dark=roman_dark_WFI01_IMAGE_STRESS_TEST_16_MA_TABLE_998_D1.asdf",
        "roman_elp",
        rtdata.input,
    ]
    ExposurePipeline.from_cmdline(args)

    # Perform DMS tests
    # Initial prep
    pipeline = ExposurePipeline()
    # rtdata.get_data(f"WFI/image/{output}")
    model = rdm.open(output, lazy_load=False)
    uncal_data = rdm.open(input_data, lazy_load=False)

    # DMS280 result is an ImageModel
    pipeline.log.info(
        "DMS414 MSG: Testing that result is a Level 2 model......."
        + passfail(isinstance(model, rdm.datamodels.ImageModel))
    )

    pipeline.log.info(
        "DMS414 MSG: Testing that there are 16 resultants in the input file......."
        + passfail(len(uncal_data.meta.exposure.read_pattern) == 16)
    )

    pipeline.log.info(
        "DMS414 MSG: Testing that there are 16 resultants listed in the output file......."
        + passfail(len(model.meta.exposure.read_pattern) == 16)
    )

    # Ensure step completion is as expected
    assert model.meta.cal_step.dq_init == "COMPLETE"
    assert model.meta.cal_step.saturation == "COMPLETE"
    assert model.meta.cal_step.linearity == "COMPLETE"
    assert model.meta.cal_step.dark == "COMPLETE"
    assert model.meta.cal_step.jump == "COMPLETE"
    assert model.meta.cal_step.ramp_fit == "COMPLETE"
    assert model.meta.cal_step.assign_wcs == "COMPLETE"
    assert model.meta.cal_step.flat_field == "SKIPPED"
    assert model.meta.cal_step.photom == "SKIPPED"
