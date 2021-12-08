import os
import pytest

from romancal.stpipe import RomanStep, RomanPipeline
from romancal.pipeline.exposure_pipeline import ExposurePipeline
import roman_datamodels as rdm

from romancal.assign_wcs.assign_wcs_step import AssignWcsStep

from .regtestdata import compare_asdf
import copy

@pytest.mark.bigdata
def test_level2_image_processing_pipeline(rtdata, ignore_asdf_paths):
    rtdata.get_data("WFI/image/l1_0001.asdf")
    rtdata.input = "l1_0001.asdf"

    # Test Pipeline
    output = "l1_0001_cal.asdf"
    rtdata.output = output
    args = ["--disable-crds-steppars",
            "--steps.jump.rejection_threshold=180.0",
            "--steps.jump.three_group_rejection_threshold=185.0",
            "roman_elp", rtdata.input,
            ]
    ExposurePipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)


    # Perform DMS tests
    # Initial prep
    model = rdm.open(rtdata.output)
    pipeline = ExposurePipeline()

    # DMS86 instrument artifact correction tests
    pipeline.log.info('Status of the step:             assign_wcs    ' +
                  str(model.meta.cal_step.assign_wcs))
    pipeline.log.info('DMS86 MSG: Testing completion of wcs assignment in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.assign_wcs == 'COMPLETE'))
    assert model.meta.cal_step.assign_wcs == 'COMPLETE'
    pipeline.log.info('Status of the step:             flat_field    ' +
                  str(model.meta.cal_step.flat_field))
    pipeline.log.info('DMS86 MSG: Testing completion of flat fielding in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.flat_field == 'COMPLETE'))
    assert model.meta.cal_step.flat_field == 'COMPLETE'
    pipeline.log.info('Status of the step:             dark          ' +
                  str(model.meta.cal_step.dark))
    pipeline.log.info('DMS86 MSG: Testing completion of dark correction in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.dark == 'COMPLETE'))
    assert model.meta.cal_step.dark == 'COMPLETE'
    pipeline.log.info('Status of the step:             dq_init       ' +
                  str(model.meta.cal_step.dq_init))
    pipeline.log.info('DMS86 MSG: Testing completion of data quality correction in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.dq_init == 'COMPLETE'))
    assert model.meta.cal_step.dq_init == 'COMPLETE'
    pipeline.log.info('Status of the step:             jump          ' +
                  str(model.meta.cal_step.jump))
    pipeline.log.info('DMS86 MSG: Testing completion of jump detection in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.jump == 'COMPLETE'))
    assert model.meta.cal_step.jump == 'COMPLETE'
    pipeline.log.info('Status of the step:             linearity     ' +
                  str(model.meta.cal_step.assign_wcs))
    pipeline.log.info('DMS86 MSG: Testing completion of linearity correction in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.linearity == 'COMPLETE'))
    assert model.meta.cal_step.linearity == 'COMPLETE'
    pipeline.log.info('Status of the step:             ramp_fit      ' +
                  str(model.meta.cal_step.ramp_fit))
    pipeline.log.info('DMS86 MSG: Testing completion of ramp fitting in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.ramp_fit == 'COMPLETE'))
    assert model.meta.cal_step.ramp_fit == 'COMPLETE'
    pipeline.log.info('Status of the step:             saturation    ' +
                  str(model.meta.cal_step.saturation))
    pipeline.log.info('DMS86 MSG: Testing completion of saturation detection in'
                  'Level 2 image output.......' +
                  str(model.meta.cal_step.saturation == 'COMPLETE'))
    assert model.meta.cal_step.saturation == 'COMPLETE'

    # DMS87 data quality tests
    pipeline.log.info('DMS87 MSG: Testing existence of data quality array (dq) in'
                  'Level 2 image output.......' +
                  str("dq" in model.keys()))
    assert "dq" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of general error array (err) in'
                  'Level 2 image output.......' +
                  str("err" in model.keys()))
    assert "err" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of Poisson noise variance array (var_poisson) in'
                  'Level 2 image output.......' +
                  str("var_poisson" in model.keys()))
    assert "var_poisson" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of read noise variance array (var_rnoise) in'
                  'Level 2 image output.......' +
                  str("var_rnoise" in model.keys()))
    assert "var_rnoise" in model.keys()
    pipeline.log.info('DMS87 MSG: Testing existence of flatfield uncertainty variance array (var_flat) in'
                  'Level 2 image output.......' +
                  str("var_flat" in model.keys()))
    assert "var_flat" in model.keys()

    # DMS88 total exposure time test
    pipeline.log.info('DMS88 MSG: Testing existence of total exposure time (exposure_time)'
                  'in Level 2 image output.......' +
                  str("exposure_time" in model.meta.exposure))
    assert "exposure_time" in model.meta.exposure

    # DMS89 WCS tests
    pipeline.log.info('DMS89 MSG: Testing that the wcs bounding'
                  'box was generated.......' +
                  str(((len(model.meta.wcs.bounding_box) == 2) and
                        (type(model.meta.wcs.bounding_box[0][0]) == float))))
    assert ((len(model.meta.wcs.bounding_box) == 2) and
            (type(model.meta.wcs.bounding_box[0][0]) == float))

    # Save original wcs information
    orig_wcs = copy.deepcopy(model.meta.wcs)
    del model.meta['wcs']

    # Create new pointing for the model
    # RA & Dec are each shifted + 10 degrees, unless they are near the upper limit,
    # in which case they are shifted -10 degrees.
    if model.meta.wcsinfo.ra_ref < 350.0:
        model.meta.wcsinfo.ra_ref += 10.0
    else:
        model.meta.wcsinfo.ra_ref -= 10.0

    if model.meta.wcsinfo.dec_ref < 80.0:
        model.meta.wcsinfo.dec_ref += 10.0
    else:
        model.meta.wcsinfo.dec_ref -= 10.0

    # Create new wcs object for the new pointing
    model = AssignWcsStep.call(model)

    pipeline.log.info('DMS89 MSG: Testing that the different pointings create differing wcs.......' +
                  str(model.meta.wcs.footprint != orig_wcs.footprint))
    assert model.meta.wcs.footprint != orig_wcs.footprint
