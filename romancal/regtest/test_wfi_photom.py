"""Regression tests for the photom step of the Roman pipeline"""
import math
import pytest

import roman_datamodels as rdm
from astropy import units as u
from romancal.stpipe import RomanStep
from romancal.step import PhotomStep
from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_absolute_photometric_calibration(rtdata, ignore_asdf_paths):
    """DMS140 Test: Testing application of photometric correction using
       CRDS selected photom file."""

    input_data = "r0000201001001001001_01101_0001_WFI01_flat.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Define step (for running and log access)
    step = PhotomStep()

    #  In Wide Field Imaging mode, the DMS shall generate Level 2 science
    # data products with absolute photometry calibrated in the WFI filter
    # used for the exposure.
    step.log.info('DMS140 MSG: Testing absolute photometric '
                  'calibrated image data. '
                  'Success is creation of a Level 2 image file with '
                  'CRDS selected photom file applied.')

    step.log.info('DMS140 MSG: Image data file: '
                  f'{rtdata.input.rsplit("/", 1)[1]}')

    # Note: if any of the following tests fail, check for a different
    # photom match from CRDS. Values come from roman_wfi_photom_0034.asdf

    # Test PhotomStep
    output = "r0000201001001001001_01101_0001_WFI01_photom.asdf"
    rtdata.output = output
    args = ["romancal.step.PhotomStep", rtdata.input]
    step.log.info('DMS140 MSG: Running photometric conversion step.'
                  ' The first ERROR is expected, due to extra CRDS parameters'
                  ' not having been implemented yet.')
    RomanStep.from_cmdline(args)
    photom_out = rdm.open(rtdata.output)

    step.log.info(f'DMS140 MSG: Photom step recorded as complete? : '
                  f'{photom_out.meta.cal_step.photom == "COMPLETE"}')
    assert photom_out.meta.cal_step.photom == "COMPLETE"

    step.log.info('DMS140 MSG: Photom megajansky conversion calculated? : ' +
                  str((photom_out.meta.photometry.conversion_megajanskys.unit == u.MJy / u.sr) and
                       math.isclose(photom_out.meta.photometry.conversion_megajanskys.value,
                                    0.3324, abs_tol=0.0001)))
    assert photom_out.meta.photometry.conversion_megajanskys.unit == u.MJy / u.sr
    assert math.isclose(photom_out.meta.photometry.conversion_megajanskys.value,
                        0.3324, abs_tol=0.0001)

    step.log.info('DMS140 MSG: Photom microjanskys conversion calculated? : ' +
                  str((photom_out.meta.photometry.conversion_microjanskys.unit ==
                       u.uJy / u.arcsec ** 2) and
                      (math.isclose(photom_out.meta.photometry.conversion_microjanskys.value,
                                    7.81320, abs_tol=0.0001))))
    assert photom_out.meta.photometry.conversion_microjanskys.unit == u.uJy / u.arcsec ** 2
    assert math.isclose(photom_out.meta.photometry.conversion_microjanskys.value,
                        7.81320, abs_tol=0.0001)

    step.log.info('DMS140 MSG: Pixel area in steradians calculated? : ' +
                  str((photom_out.meta.photometry.pixelarea_steradians.unit == u.sr) and
                      (math.isclose(photom_out.meta.photometry.pixelarea_steradians.value,
                                    2.8083e-13, abs_tol=1.0e-17))))
    assert photom_out.meta.photometry.pixelarea_steradians.unit == u.sr
    assert math.isclose(photom_out.meta.photometry.pixelarea_steradians.value, 2.8083e-13,
                        abs_tol=1.0e-17)

    step.log.info('DMS140 MSG: Pixel area in square arcseconds calculated? : ' +
                  str((photom_out.meta.photometry.pixelarea_arcsecsq.unit == u.arcsec ** 2) and
                      (math.isclose(photom_out.meta.photometry.pixelarea_arcsecsq.value,
                                    0.011948, abs_tol=1.0e-6))))
    assert photom_out.meta.photometry.pixelarea_arcsecsq.unit == u.arcsec ** 2
    assert math.isclose(photom_out.meta.photometry.pixelarea_arcsecsq.value, 0.011948,
                        abs_tol=1.0e-6)

    step.log.info('DMS140 MSG: Photom megajansky conversion uncertainty calculated? : ' +
                  str((photom_out.meta.photometry.conversion_megajanskys_uncertainty.unit
                       == u.MJy / u.sr) and
                      (math.isclose(photom_out.meta.photometry.\
                                    conversion_megajanskys_uncertainty.value, 0.0,
                                    abs_tol=1.0e-6))))
    assert photom_out.meta.photometry.conversion_megajanskys_uncertainty.unit == u.MJy / u.sr
    assert math.isclose(photom_out.meta.photometry.conversion_megajanskys_uncertainty.value, 0.0,
                        abs_tol=1.0e-6)

    step.log.info('DMS140 MSG: Photom megajansky conversion uncertainty calculated? : ' +
                  str((photom_out.meta.photometry.conversion_microjanskys_uncertainty.unit ==
                       u.uJy / u.arcsec ** 2) and
                      (math.isclose(photom_out.meta.photometry.\
                                    conversion_microjanskys_uncertainty.value, 0.0,
                                    abs_tol=1.0e-6))))
    assert photom_out.meta.photometry.conversion_microjanskys_uncertainty.unit ==  \
           u.uJy / u.arcsec ** 2
    assert math.isclose(photom_out.meta.photometry.conversion_microjanskys_uncertainty.value, 0.0,
                        abs_tol=1.0e-6)

    rtdata.get_truth(f"truth/WFI/image/{output}")
    step.log.info(f'DMS140 MSG: Was the proper absolute photometry calibrated image data produced?'
                  f' : {(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}')
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None
