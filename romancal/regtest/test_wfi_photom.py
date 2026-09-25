"""Regression tests for the photom step of the Roman pipeline"""

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

    # assign_wcs stage is input to photom
    input_data = "r0000101001001001001_0001_wfi01_f158_assignwcs.asdf"
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

    # Test PhotomStep
    output = "r0000101001001001001_0001_wfi01_f158_photom.asdf"
    rtdata.output = output
    # Fetch the truth before any assertions so that okify files can
    # be generated
    rtdata.get_truth(f"truth/WFI/image/{output}")
    args = ["romancal.step.PhotomStep", rtdata.input]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    with rdm.open(rtdata.output) as photom_out:
        dms_logger.info(
            "DMS140 MSG: Photom step recorded as complete? :"
            f" {photom_out.meta.cal_step.photom == 'COMPLETE'}"
        )
        assert photom_out.meta.cal_step.photom == "COMPLETE"

        # check for reasonable values of conversion_megajansky,
        # pixel_area, and uncertainty
        photometry = photom_out.meta.photometry
        conv = photometry.conversion_megajanskys
        conv_ok = 0.5 < conv < 1.0
        dms_logger.info(
            f"DMS140 MSG: Photom megajansky conversion calculated? : {conv_ok}"
        )
        assert conv_ok, f"conversion_megajanskys = {conv}"

        # nominal WFI pixel is 0.11 arcsec, or ~2.84e-13 sr
        area = photometry.pixel_area
        area_ok = 2.5e-13 < area < 3.2e-13
        dms_logger.info(f"DMS140 MSG: Pixel area in steradians calculated? : {area_ok}")
        assert area_ok, f"pixel_area = {area}"

        unc = photometry.conversion_megajanskys_uncertainty
        unc_ok = 0 < unc < 0.2 * conv
        dms_logger.info(
            f"DMS140 MSG: Photom megajansky conversion uncertainty calculated? : {unc_ok}"
        )
        assert unc_ok, f"conversion_megajanskys_uncertainty = {unc}"

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    dms_logger.info(
        "DMS140 MSG: Was the proper absolute photometry calibrated image data produced?"
        f" : {diff.identical}"
    )
    assert diff.identical, diff.report()
