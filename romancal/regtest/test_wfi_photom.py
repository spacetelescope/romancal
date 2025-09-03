"""Regression tests for the photom step of the Roman pipeline"""

import math

import pytest
import roman_datamodels as rdm

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_absolute_photometric_calibration(
    rtdata, ignore_asdf_paths, resource_tracker, request, dms_logger
):
    """DMS140 Test: Testing application of photometric correction using
    CRDS selected photom file."""

    input_data = "r0000101001001001001_0001_wfi01_f158_flat.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    #  In Wide Field Imaging mode, the DMS shall generate Level 2 science
    # data products with absolute photometry calibrated in the WFI filter
    # used for the exposure.
    dms_logger.info(
        "DMS140 MSG: Testing absolute photometric "
        "calibrated image data. "
        "Success is creation of a Level 2 image file with "
        "CRDS selected photom file applied."
    )

    dms_logger.info(f"DMS140 MSG: Image data file: {rtdata.input.rsplit('/', 1)[1]}")

    # Note: if any of the following tests fail, check for a different
    # photom match from CRDS. Values come from roman_wfi_photom_0034.asdf

    # Test PhotomStep
    output = "r0000101001001001001_0001_wfi01_f158_photom.asdf"
    rtdata.output = output
    args = ["romancal.step.PhotomStep", rtdata.input]
    dms_logger.info(
        "DMS140 MSG: Running photometric conversion step."
        " The first ERROR is expected, due to extra CRDS parameters"
        " not having been implemented yet."
    )
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    photom_out = rdm.open(rtdata.output)

    dms_logger.info(
        "DMS140 MSG: Photom step recorded as complete? :"
        f" {photom_out.meta.cal_step.photom == 'COMPLETE'}"
    )
    assert photom_out.meta.cal_step.photom == "COMPLETE"

    convval = 0.73678
    dms_logger.info(
        "DMS140 MSG: Photom megajansky conversion calculated? : "
        + str(
            math.isclose(
                photom_out.meta.photometry.conversion_megajanskys,
                convval,
                abs_tol=0.0001,
            )
        )
    )
    assert math.isclose(
        photom_out.meta.photometry.conversion_megajanskys, convval, abs_tol=0.0001
    )

    dms_logger.info(
        "DMS140 MSG: Pixel area in steradians calculated? : "
        + str(
            math.isclose(
                photom_out.meta.photometry.pixel_area,
                2.8083e-13,
                abs_tol=1.0e-17,
            )
        )
    )
    assert math.isclose(
        photom_out.meta.photometry.pixel_area,
        2.8083e-13,
        abs_tol=1.0e-17,
    )

    uncval = 0.02866405
    dms_logger.info(
        "DMS140 MSG: Photom megajansky conversion uncertainty calculated? : "
        + str(
            math.isclose(
                photom_out.meta.photometry.conversion_megajanskys_uncertainty,
                uncval,
                abs_tol=1.0e-6,
            )
        )
    )
    assert math.isclose(
        photom_out.meta.photometry.conversion_megajanskys_uncertainty,
        uncval,
        abs_tol=1.0e-6,
    )

    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    dms_logger.info(
        "DMS140 MSG: Was the proper absolute photometry calibrated image data produced?"
        f" : {diff.identical}"
    )
    assert diff.identical, diff.report()
