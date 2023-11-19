""" Module to test rampfit with optional output """
from metrics_logger.decorators import metrics_logger
from pathlib import Path
import pytest

import roman_datamodels as rdm
from romancal.lib.dms import log_result
from romancal.lib.suffix import replace_suffix
from romancal.ramp_fitting import RampFitStep
from romancal.regtest.conftest import ignore_asdf_paths, rtdata_module
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


# ##########
# Conditions
# ##########
def cond_is_asdf(requirement, model, expected_path):
    """Check that the filename has the correct file type"""
    result = expected_path.exists() and expected_path.suffix == '.asdf'
    log_result(requirement, 'Testing that result file path has file type "asdf"', result)
    return result


def cond_is_imagemodel(requirement, model, expected_path):
    """Check that the result is an ImageModel"""
    result = isinstance(model, rdm.datamodels.ImageModel)
    log_result(requirement, 'Testing that the result model is Level 2', result)
    return result


def cond_is_rampfit(requirement, model, expected_path):
    """Check that the calibration suffix is 'rampfit'"""
    result = expected_path.exists() and expected_path.stem.endswith('rampfit')
    log_result(requirement, 'Testing that the result file has the suffix "rampfit"', result)
    return result


def cond_is_step_complete(requirement, model, expected_path):
    """Check that the calibration step is marked complete"""
    result = model.meta.cal_step.ramp_fit == 'COMPLETE'
    log_result(requirement, 'Testing that RampFitStep completed', result)
    return result


def cond_is_uneven(requirement, model, expected_path):
    """Verify that the provided model represents uneven ramps

    Parameters
    ----------
    rampfit_result : `roman_datamodels.ImageModel`
        Model created from `RampFitStep`
    """
    length_set = {len(resultant) for resultant in model.meta.exposure.read_pattern}

    result = len(length_set) > 1
    log_result(requirement, 'Testing that the ramps are uneven', result)
    return result


def cond_science_verification(requirement, model, expected_path, rtdata_module, ignore_asdf_paths):
    """Check against expected data results"""
    diff = compare_asdf(rtdata_module.output, rtdata_module.truth, **ignore_asdf_paths)

    result = diff.identical
    if not result:
        diff.report()
    log_result(requirement, 'Testing science veracity', result)
    return result


CONDITIONS_FULL = [cond_is_asdf, cond_is_imagemodel, cond_is_rampfit, cond_is_step_complete, cond_is_uneven]

# ######################
# fixtures and utilities
# ######################
@pytest.fixture(scope='module',
                params=[('DMS362', Path('WFI/image/r0000101001001001001_01101_0001_WFI01_dqinit.asdf'), CONDITIONS_FULL),
                        ('DMS366', Path('WFI/grism/r0000201001001001001_01101_0001_WFI05_dqinit.asdf'), CONDITIONS_FULL),
                        ('DMS363', Path('WFI/image/r0000101001001001001_01101_0003_WFI01_dqinit.asdf'), CONDITIONS_FULL),
                        ('DMS367', Path('WFI/grism/r0000201001001001001_01101_0003_WFI05_dqinit.asdf'), CONDITIONS_FULL)])
def rampfit_result(request, rtdata_module):
    """Run RampFitStep

    Parameters
    ----------
    rtdata_module : pytest.fixture
        artifactory fixture for data retrieval

    Returns
    -------
    model, path : `ImageModel`, `pathlib.Path`
        Model and path to model.
    """
    # Setup inputs
    requirement, artifactory_path, conditions = request.param
    input_data = rtdata_module.get_data(str(artifactory_path))
    rtdata_module.input = input_data

    # Execute the step
    result_model = RampFitStep.call(input_data, save_results=True)

    # Setup outputs
    input_data_path = Path(input_data)
    output = (replace_suffix(input_data_path.stem, 'rampfit') + '.asdf')
    expected_path = input_data_path.parent / output

    # Get truths
    rtdata_module.output = output
    output_artifactory_path = Path('truth') / artifactory_path.parent / output
    rtdata_module.get_truth(str(output_artifactory_path))

    try:
        yield requirement, result_model, expected_path, conditions
    finally:
        result_model.close()


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


# #####
# Tests
# #####
@pytest.mark.bigdata
def test_rampfit_success(rampfit_result, rtdata_module, ignore_asdf_paths):
    """Test rampfit result against various conditions"""
    requirement, model, expected_path, conditions = rampfit_result
    success = True
    for condition in conditions:
        success = success and condition(requirement, model, expected_path)

    # Always do a full regression check.
    success = success and cond_science_verification(requirement, model, expected_path, rtdata_module, ignore_asdf_paths)

    @metrics_logger(requirement)
    def test_success():
        assert success
    test_success()



@pytest.mark.bigdata
def test_ramp_fitting_step(rtdata, ignore_asdf_paths):
    """Testing the ramp fitting step"""
    input_data = "r0000101001001001001_01101_0001_WFI01_darkcurrent.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    args = [
        "romancal.step.RampFitStep",
        rtdata.input,
    ]
    RomanStep.from_cmdline(args)
    output = "r0000101001001001001_01101_0001_WFI01_rampfit.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()
