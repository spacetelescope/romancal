""" Roman tests for flat field correction """
import copy
import pytest

import numpy as np
from numpy.testing import assert_allclose

import roman_datamodels as rdm
from romancal.pipeline.exposure_pipeline import ExposurePipeline
from romancal.assign_wcs.assign_wcs_step import AssignWcsStep
from .regtestdata import compare_asdf

from gwcs.wcstools import grid_from_bounding_box


def passfail(bool_expr):
    """ set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


@pytest.mark.bigdata
@pytest.mark.soctests
def test_level2_image_processing_pipeline(rtdata, ignore_asdf_paths):
    """Tests for flat field imaging processing requirements DMS86 & DMS 87 """
    input_data = "r0000101001001001001_01101_0001_WFI01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000101001001001001_01101_0001_WFI01_cal.asdf"
    rtdata.output = output
    args = ["--disable-crds-steppars",
            "--steps.jump.rejection_threshold=180.0",
            "--steps.jump.three_group_rejection_threshold=185.0",
            "--steps.jump.four_group_rejection_threshold=190",
            "roman_elp", rtdata.input,
            ]
    ExposurePipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None

    # Perform DMS tests
    # Initial prep
    model = rdm.open(rtdata.output)
    pipeline = ExposurePipeline()

    # DMS86 instrument artifact correction tests
    pipeline.log.info('Status of the step:             assign_wcs    ' +
                      str(model.meta.cal_step.assign_wcs))
    pipeline.log.info('DMS86 MSG: Testing completion of wcs assignment in'
                      'Level 2 image output.......' +
                      passfail(model.meta.cal_step.assign_wcs == 'COMPLETE'))
    assert model.meta.cal_step.assign_wcs == 'COMPLETE'
    pipeline.log.info('Status of the step:             flat_field    ' +
                      str(model.meta.cal_step.flat_field))
    pipeline.log.info('DMS86 MSG: Testing completion of flat fielding in'
                      'Level 2 image output.......' +
                      passfail(model.meta.cal_step.flat_field == 'PASS'))
    assert model.meta.cal_step.flat_field == 'COMPLETE'
    pipeline.log.info('Status of the step:             dark          ' +
                      str(model.meta.cal_step.dark))
    pipeline.log.info('DMS86 MSG: Testing completion of dark correction in'
                      'Level 2 image output.......' +
                      passfail(model.meta.cal_step.dark == 'COMPLETE'))
    assert model.meta.cal_step.dark == 'COMPLETE'
    pipeline.log.info('Status of the step:             dq_init       ' +
                      str(model.meta.cal_step.dq_init))
    pipeline.log.info('DMS86 MSG: Testing completion of data quality'
                      ' correction in Level 2 image output.......' +
                      passfail(model.meta.cal_step.dq_init == 'COMPLETE'))
    assert model.meta.cal_step.dq_init == 'COMPLETE'
    pipeline.log.info('Status of the step:             jump          ' +
                      str(model.meta.cal_step.jump))
    pipeline.log.info('DMS86 MSG: Testing completion of jump detection in'
                      'Level 2 image output.......' +
                      passfail(model.meta.cal_step.jump == 'COMPLETE'))
    assert model.meta.cal_step.jump == 'COMPLETE'
    pipeline.log.info('Status of the step:             linearity     ' +
                      str(model.meta.cal_step.linearity))
    pipeline.log.info('DMS86 MSG: Testing completion of linearity correction'
                      ' in Level 2 image output.......' +
                      passfail(model.meta.cal_step.linearity == 'COMPLETE'))
    assert model.meta.cal_step.linearity == 'COMPLETE'
    pipeline.log.info('Status of the step:             ramp_fit      ' +
                      str(model.meta.cal_step.ramp_fit))
    pipeline.log.info('DMS86 MSG: Testing completion of ramp fitting in'
                      'Level 2 image output.......' +
                      passfail(model.meta.cal_step.ramp_fit == 'COMPLETE'))
    assert model.meta.cal_step.ramp_fit == 'COMPLETE'
    pipeline.log.info('Status of the step:             saturation    ' +
                      str(model.meta.cal_step.saturation))
    pipeline.log.info('DMS86 MSG: Testing completion of saturation detection '
                      'in Level 2 image output.......' +
                      passfail(model.meta.cal_step.saturation == 'COMPLETE'))
    assert model.meta.cal_step.saturation == 'COMPLETE'

    # DMS-129 tests
    # check if assign_wcs step is complete
    pipeline.log.info('DMS-129 MSG: Status of the step:             assign_wcs    ' +
                    str(model.meta.cal_step.assign_wcs))
    pipeline.log.info('DMS-129 MSG: Testing completion of WCS assignment in'
                    'Level 2 image output.......' +
                    passfail(model.meta.cal_step.assign_wcs == 'COMPLETE'))
    assert model.meta.cal_step.assign_wcs == 'COMPLETE'
    # check if WCS exists
    pipeline.log.info('DMS-129 MSG: Testing that a WCS object exists    ')
    pipeline.log.info('DMS-129 MSG: Testing that WCS exists in'
                    'Level 2 image output.......' +
                    passfail(model.meta.wcs is not None))
    assert model.meta.wcs is not None
    pipeline.log.info('DMS-129 MSG: Testing that geometric distortion information is available in'
                    'Level 2 image output.......' +
                    passfail('v2v3' in model.meta.wcs.available_frames))
    assert 'v2v3' in model.meta.wcs.available_frames
    # compare coordinates before and after distortion correction has been applied
    # 1 - get new image array based on the model
    x0, y0 = grid_from_bounding_box(model.meta.wcs.bounding_box)
    # 2 - apply the distortion-corrected WCS solution to new image array
    corrected_coords = model.meta.wcs(x0, y0)
    # 3 - apply the transformation from 'v2v3' to 'world' without distortion correction
    original_coords = model.meta.wcs.get_transform('v2v3', 'world')(x0, y0)
    # compare both results to make sure they don't match
    # (which means the distortion correction was actually applied to the model)
    pipeline.log.info('DMS-129 MSG: Testing that distortion correction was applied to'
                'Level 2 image output.......' +
                passfail((corrected_coords[0] != original_coords[0]).all() & (corrected_coords[1] != original_coords[1]).all()))
    assert (corrected_coords[0] != original_coords[0]).all()
    assert (corrected_coords[1] != original_coords[1]).all()

    # DMS87 data quality tests
    pipeline.log.info('DMS87 MSG: Testing existence of data quality array (dq)'
                      ' in Level 2 image output.......' +
                      passfail("dq" in model.keys()))
    assert "dq" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of general error array '
                      '(err) in Level 2 image output.......' +
                      passfail("err" in model.keys()))
    assert "err" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of Poisson noise variance'
                      'array (var_poisson) in Level 2 image output.......' +
                      passfail("var_poisson" in model.keys()))
    assert "var_poisson" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of read noise variance '
                      'array (var_rnoise) in level 2 image output.......' +
                      passfail("var_rnoise" in model.keys()))
    assert "var_rnoise" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of flatfield uncertainty '
                      'variance array (var_flat) in Level 2 image output....' +
                      passfail("var_flat" in model.keys()))
    assert "var_flat" in model.keys()

    # DMS88 total exposure time test
    pipeline.log.info('DMS88 MSG: Testing existence of total exposure time '
                      '(exposure_time) in Level 2 image output.......' +
                      passfail("exposure_time" in model.meta.exposure))
    assert "exposure_time" in model.meta.exposure

    # DMS89 WCS tests
    pipeline.log.info('DMS89 MSG: Testing that the wcs bounding'
                      'box was generated.......' +
                      passfail((len(model.meta.wcs.bounding_box) == 2)))
    assert len(model.meta.wcs.bounding_box) == 2

    # Save original wcs information
    orig_wcs = copy.deepcopy(model.meta.wcs)
    del model.meta['wcs']

    # Create new pointing for the model
    # RA & Dec are each shifted + 10 degrees, unless they are near
    # the upper limit, in which case they are shifted -10 degrees.
    delta = [10.0, 10.0]
    if model.meta.wcsinfo.ra_ref >= 350.0:
        delta[0] *= -1.0

    if model.meta.wcsinfo.dec_ref >= 80.0:
        delta[1] *= -1.0

    model.meta.wcsinfo.ra_ref += delta[0]
    model.meta.wcsinfo.dec_ref += delta[1]

    # Create new wcs object for the new pointing
    model = AssignWcsStep.call(model)

    rtdata.output = output.rsplit(".", 1)[0] + "_repoint.asdf"
    model.to_asdf(rtdata.output)

    # Test that repointed file matches truth
    rtdata.get_truth("truth/WFI/image/" + output.rsplit(".", 1)[0] +
                     "_repoint.asdf")
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None

    pipeline.log.info('DMS89 MSG: Testing that the different pointings '
                      'create differing wcs.......'
                      + passfail(((np.abs(orig_wcs(2048, 2048)[0] -
                                          model.meta.wcs(2048, 2048)[0])) -
                                  10.0) < 1.0))
    assert_allclose([orig_wcs(2048, 2048)[0] +
                    delta[0], orig_wcs(2048, 2048)[1] + delta[1]],
                    model.meta.wcs(2048, 2048), atol=1.0)


@pytest.mark.bigdata
@pytest.mark.soctests
def test_level2_grism_processing_pipeline(rtdata, ignore_asdf_paths):
    """Tests for flat field grism processing requirements DMS90 & DMS 91 """
    input_data = "r0000201001001001002_01101_0001_WFI01_uncal.asdf"
    rtdata.get_data(f"WFI/grism/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000201001001001002_01101_0001_WFI01_cal.asdf"
    rtdata.output = output
    args = ["--disable-crds-steppars",
            "--steps.jump.rejection_threshold=180.0",
            "--steps.jump.three_group_rejection_threshold=185.0",
            "--steps.jump.four_group_rejection_threshold=190",
            "roman_elp", rtdata.input,
            ]
    ExposurePipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/grism/{output}")
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None

    # Perform DMS tests
    # Initial prep
    model = rdm.open(rtdata.output)
    pipeline = ExposurePipeline()

    # DMS90 instrument artifact correction tests
    pipeline.log.info('Status of the step:             assign_wcs    ' +
                      str(model.meta.cal_step.assign_wcs))
    pipeline.log.info('DMS90 MSG: Testing completion of wcs assignment in'
                      'Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.assign_wcs == 'COMPLETE'))
    assert model.meta.cal_step.assign_wcs == 'COMPLETE'
    pipeline.log.info('Status of the step:             flat_field    ' +
                      str(model.meta.cal_step.flat_field))
    pipeline.log.info('DMS90 MSG: Testing expected skip of flat fielding in'
                      'Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.flat_field == 'SKIPPED'))
    assert model.meta.cal_step.flat_field == 'SKIPPED'
    pipeline.log.info('Status of the step:             dark          ' +
                      str(model.meta.cal_step.dark))
    pipeline.log.info('DMS90 MSG: Testing completion of dark correction in'
                      'Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.dark == 'COMPLETE'))
    assert model.meta.cal_step.dark == 'COMPLETE'
    pipeline.log.info('Status of the step:             dq_init       ' +
                      str(model.meta.cal_step.dq_init))
    pipeline.log.info('DMS90 MSG: Testing completion of data quality '
                      'correction in Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.dq_init == 'COMPLETE'))
    assert model.meta.cal_step.dq_init == 'COMPLETE'
    pipeline.log.info('Status of the step:             jump          ' +
                      str(model.meta.cal_step.jump))
    pipeline.log.info('DMS90 MSG: Testing completion of jump detection in'
                      'Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.jump == 'COMPLETE'))
    assert model.meta.cal_step.jump == 'COMPLETE'
    pipeline.log.info('Status of the step:             linearity     ' +
                      str(model.meta.cal_step.assign_wcs))
    pipeline.log.info('DMS90 MSG: Testing completion of linearity correction '
                      'in Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.linearity == 'COMPLETE'))
    assert model.meta.cal_step.linearity == 'COMPLETE'
    pipeline.log.info('Status of the step:             ramp_fit      ' +
                      str(model.meta.cal_step.ramp_fit))
    pipeline.log.info('DMS90 MSG: Testing completion of ramp fitting in'
                      'Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.ramp_fit == 'COMPLETE'))
    assert model.meta.cal_step.ramp_fit == 'COMPLETE'
    pipeline.log.info('Status of the step:             saturation    ' +
                      str(model.meta.cal_step.saturation))
    pipeline.log.info('DMS90 MSG: Testing completion of saturation detection '
                      'in Level 2 spectroscopic output.......' +
                      passfail(model.meta.cal_step.saturation == 'COMPLETE'))
    assert model.meta.cal_step.saturation == 'COMPLETE'

    # DMS91 data quality tests
    pipeline.log.info('DMS91 MSG: Testing existence of data quality array (dq)'
                      ' in Level 2 spectroscopic output.......' +
                      passfail("dq" in model.keys()))
    assert "dq" in model.keys()
    pipeline.log.info('DMS91 MSG: Testing existence of general error '
                      'array (err) in Level 2 spectroscopic output.......' +
                      passfail("err" in model.keys()))
    assert "err" in model.keys()
    pipeline.log.info('DMS91 MSG: Testing existence of Poisson noise variance '
                      'array (var_poisson) in Level 2 spectroscopic output..' +
                      passfail("var_poisson" in model.keys()))
    assert "var_poisson" in model.keys()
    pipeline.log.info('DMS91 MSG: Testing existence of read noise variance '
                      'array (var_rnoise) in Level 2 spectroscopic output...' +
                      passfail("var_rnoise" in model.keys()))
    assert "var_rnoise" in model.keys()

    # DMS88 total exposure time test
    pipeline.log.info('DMS88 MSG: Testing existence of total exposure '
                      'time (exposure_time) in Level 2 spectroscopic output.' +
                      passfail("exposure_time" in model.meta.exposure))
    assert "exposure_time" in model.meta.exposure

    # DMS93 WCS tests
    pipeline.log.info('DMS93 MSG: Testing that the wcs bounding'
                      'box was generated.......' +
                      passfail((len(model.meta.wcs.bounding_box) == 2)))
    assert len(model.meta.wcs.bounding_box) == 2

    # Save original wcs information
    orig_wcs = copy.deepcopy(model.meta.wcs)
    del model.meta['wcs']

    # Create new pointing for the model
    # RA & Dec are each shifted + 10 degrees, unless they are near
    # the upper limit, in which case they are shifted -10 degrees.
    delta = [10.0, 10.0]
    if model.meta.wcsinfo.ra_ref >= 350.0:
        delta[0] *= -1.0

    if model.meta.wcsinfo.dec_ref >= 80.0:
        delta[1] *= -1.0

    model.meta.wcsinfo.ra_ref += delta[0]
    model.meta.wcsinfo.dec_ref += delta[1]

    # Create new wcs object for the new pointing
    model = AssignWcsStep.call(model)

    rtdata.output = output.rsplit(".", 1)[0] + "_repoint.asdf"
    model.to_asdf(rtdata.output)

    # Test that repointed file matches truth
    rtdata.get_truth("truth/WFI/grism/" + output.rsplit(".", 1)[0] +
                     "_repoint.asdf")
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None

    pipeline.log.info('DMS93 MSG: Testing that the different pointings '
                      'create differing wcs.......'
                      + passfail(((np.abs(orig_wcs(2048, 2048)[0] -
                                          model.meta.wcs(2048, 2048)[0])) -
                                  10.0) < 1.0))
    assert_allclose([orig_wcs(2048, 2048)[0] + delta[0],
                    orig_wcs(2048, 2048)[1] + delta[1]],
                    model.meta.wcs(2048, 2048),
                    atol=1.0)

@pytest.mark.bigdata

def test_processing_pipeline_all_saturated(rtdata, ignore_asdf_paths):
    """ Tests for fully saturated data skipping steps in the pipeline """
    input_data = "r0000101001001001001_01101_0001_WFI01_ALL_SATURATED_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000101001001001001_01101_0001_WFI01_ALL_SATURATED_cal.asdf"
    rtdata.output = output
    args = ["--disable-crds-steppars",
            "roman_elp", rtdata.input,
            ]
    ExposurePipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")

    assert compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None

    # Ensure step completion is as expected
    model = rdm.open(rtdata.output)

    assert model.meta.cal_step.dq_init == "COMPLETE"
    assert model.meta.cal_step.saturation == "COMPLETE"
    assert model.meta.cal_step.linearity == "SKIPPED"
    assert model.meta.cal_step.dark == "SKIPPED"
    assert model.meta.cal_step.jump == "SKIPPED"
    assert model.meta.cal_step.ramp_fit == "SKIPPED"
    assert model.meta.cal_step.assign_wcs == "SKIPPED"
    assert model.meta.cal_step.flat_field == "SKIPPED"
    assert model.meta.cal_step.photom == "SKIPPED"
