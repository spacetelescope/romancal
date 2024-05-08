""" Roman tests for the High Level Pipeline """

import os

import pytest
import roman_datamodels as rdm
from metrics_logger.decorators import metrics_logger

from romancal.pipeline.highlevel_pipeline import HighLevelPipeline

from .regtestdata import compare_asdf


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS356")
def test_level3_hlp_pipeline(rtdata, ignore_asdf_paths):
    """Tests for level 3 processing requirements DMS356"""

    cal_files = [
        "WFI/image/r0000101001001001001_01101_0001_WFI01_cal.asdf",
        "WFI/image/r0000101001001001001_01101_0002_WFI01_cal.asdf",
        "WFI/image/r0000101001001001001_01101_0003_WFI01_cal.asdf",
    ]

    for cal_file in cal_files:
        rtdata.get_data(cal_file)

    input_asn = "L3_regtest_asn.json"
    rtdata.get_data(f"WFI/image/{input_asn}")
    rtdata.input = input_asn

    # Test Pipeline
    output = "r0099101001001001001_F158_visit_i2d.asdf"
    rtdata.output = output
    args = [
        "--disable-crds-steppars",
        "roman_hlp",
        rtdata.input,
    ]
    HighLevelPipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

    # Generate thumbnail image
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    thumbnail_file = "r0099101001001001001_F158_visit_thumb.png"
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    thumbnail_file = "r0099101001001001001_F158_visit_thumb.png"

    preview_cmd = f"stpreview to {input_file} {thumbnail_file} 256 256 roman"
    os.system(preview_cmd)  # nosec

    # Generate preview image
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    preview_file = "r0099101001001001001_F158_visit_preview.png"
    preview_cmd = f"stpreview to {input_file} {preview_file} 1080 1080 roman"
    os.system(preview_cmd)  # nosec

    # Perform DMS tests
    # Initial prep
    model = rdm.open(rtdata.output, lazy_load=False)
    pipeline = HighLevelPipeline()

    # DMS356 result is an ImageModel
    pipeline.log.info(
        "DMS356 MSG: Testing that result is a Level 3 mosaic model......."
        + passfail(isinstance(model, rdm.datamodels.MosaicModel))
    )

    # DMS356 Test that skymatch step is complete
    pipeline.log.info(
        "Status of the step:             skymatch    "
        + str(model.meta.cal_step.skymatch)
    )
    # DMS356 Test that the thumbnail image exists
    pipeline.log.info(
        "Status of the step:             thumbnail image    "
        + passfail(os.path.isfile(thumbnail_file))
    )
    # DMS356 Test that the preview image exists
    pipeline.log.info(
        "Status of the step:             preview image    "
        + passfail(os.path.isfile(preview_file))
    )
    pipeline.log.info(
        "DMS86 MSG: Testing completion of skymatch in the Level 3  output......."
        + passfail(model.meta.cal_step.skymatch == "COMPLETE")
    )
    assert model.meta.cal_step.skymatch == "COMPLETE"
    pipeline.log.info(
        "Status of the step:             skymatch    "
        + str(model.meta.cal_step.skymatch)
    )
    pipeline.log.info(
        "DMS86 MSG: Testing completion of outlier detection in the Level 3 image output......."
        + passfail(model.meta.cal_step.outlier_detection == "COMPLETE")
    )
    assert model.meta.cal_step.outlier_detection == "COMPLETE"
    pipeline.log.info(
        "Status of the step:             outlier_detection    "
        + str(model.meta.cal_step.outlier_detection)
    )
    assert model.meta.cal_step.resample == "COMPLETE"
    pipeline.log.info(
        "Status of the step:             resample          "
        + str(model.meta.cal_step.resample)
    )
    pipeline.log.info(
        "DMS86 MSG: Testing completion of resample in the Level 3 image output......."
        + passfail(model.meta.cal_step.resample == "COMPLETE")
    )
