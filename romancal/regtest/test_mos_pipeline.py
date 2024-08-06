""" Roman tests for the High Level Pipeline """

import os

import pytest
import roman_datamodels as rdm
from metrics_logger.decorators import metrics_logger

from romancal.pipeline.mosaic_pipeline import MosaicPipeline

from .regtestdata import compare_asdf


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS356", "DMS374")
def test_level3_mos_pipeline(rtdata, ignore_asdf_paths):
    """Tests for level 3 processing requirements DMS356"""
    rtdata.get_asn("WFI/image/L3_regtest_asn.json")

    # Test Pipeline
    output = "r0099101001001001001_F158_visit_i2d.asdf"
    rtdata.output = output
    args = [
        "roman_mos",
        rtdata.input,
    ]
    MosaicPipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

    # Generate thumbnail image
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    thumbnail_file = "r0099101001001001001_F158_visit_thumb.png"

    preview_cmd = f"stpreview to {input_file} {thumbnail_file} 256 256 roman"
    os.system(preview_cmd)  # nosec

    # Generate preview image
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    preview_file = "r0099101001001001001_F158_visit_preview.png"
    preview_cmd = f"stpreview to {input_file} {preview_file} 1080 1080 roman"
    os.system(preview_cmd)  # nosec

    # expected catalog and segmentation files
    catalog_file = "r0099101001001001001_F158_visit_cat.asdf"
    segm_file = "r0099101001001001001_F158_visit_segm.asdf"

    # Perform DMS tests
    # Initial prep
    model = rdm.open(rtdata.output)
    pipeline = MosaicPipeline()

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
    # DMS374 Test that the output catalog exists
    pipeline.log.info(
        "Check that the catalog file exists   " + passfail(os.path.isfile(catalog_file))
    )
    # DMS374 Test that the segmentation file exists
    pipeline.log.info(
        "Check that the degmentation file exists   "
        + passfail(os.path.isfile(segm_file))
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


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS373")
def test_hlp_mosaic_pipeline(rtdata, ignore_asdf_paths):
    """Tests for level 3 mosaic requirements DMS373"""
    rtdata.get_asn("WFI/image/L3_mosaic_asn.json")

    # Test Pipeline
    output = "r0099101001001001001_r274dp63x31y81_prompt_F158_i2d.asdf"
    rtdata.output = output
    args = [
        "roman_mos",
        rtdata.input,
    ]
    MosaicPipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    pipeline = MosaicPipeline()
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

    model = rdm.open(rtdata.output, lazy_load=False)

    pipeline.log.info(
        "DMS373 MSG: Testing the creation of a Level 3 mosaic image resampled to a skycell"
        + passfail(model.meta.cal_step.resample == "COMPLETE")
    )
