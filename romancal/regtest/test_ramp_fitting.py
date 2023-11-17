""" Module to test rampfit with optional output """
from metrics_logger.decorators import metrics_logger
from pathlib import Path
import pytest

import roman_datamodels as rdm
from romancal.lib.dms import result_logger
from romancal.lib.suffix import replace_suffix
from romancal.ramp_fitting import RampFitStep
from romancal.regtest.conftest import ignore_asdf_paths, rtdata_module
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


# ######################
# fixtures and utilities
# ######################
@pytest.fixture(scope='module',
                params=[('DMS362', Path('WFI/image/r0000101001001001001_01101_0001_WFI01_dqinit.asdf')),
                        ('DMS366', Path('WFI/grism/r0000201001001001001_01101_0001_WFI05_dqinit.asdf')),
                        ('DMS363', Path('WFI/image/r0000101001001001001_01101_0003_WFI01_dqinit.asdf'))])
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
    requirement, artifactory_path = request.param
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
        yield requirement, result_model, expected_path
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
def test_is_asdf(rampfit_result):
    """Check that the filename has the correct file type"""
    requirement, _, result_path = rampfit_result

    @result_logger(requirement, 'Testing that the result file is of type "asdf"')
    def test():
        assert result_path.exists() and result_path.suffix == '.asdf'
    test()


@pytest.mark.bigdata
def test_is_imagemodel(rampfit_result):
    """Check that the result is an ImageModel"""
    requirement, model, _ = rampfit_result

    @result_logger(requirement, 'Testing that the result model is Level 2')
    def test():
        assert isinstance(model, rdm.datamodels.ImageModel)
    test()


@pytest.mark.bigdata
def test_is_rampfit(rampfit_result):
    """Check that the calibration suffix is 'rampfit'"""
    requirement, _, result_path = rampfit_result

    @result_logger(requirement, 'Testing that the result file has the suffix "rampfit"')
    def test():
        assert result_path.exists() and result_path.stem.endswith('rampfit')
    test()


@pytest.mark.bigdata
def test_is_step_complete(rampfit_result):
    """Check that the calibration step is marked complete"""
    requirement, model, _ = rampfit_result

    @result_logger(requirement, 'Testing that RampFitStep completed')
    def test():
        assert model.meta.cal_step.ramp_fit == 'COMPLETE'
    test()


@pytest.mark.bigdata
def test_is_uneven(rampfit_result):
    """Verify that the provided model represents uneven ramps

    Parameters
    ----------
    rampfit_result : `roman_datamodels.ImageModel`
        Model created from `RampFitStep`
    """
    requirement, model, _ = rampfit_result
    length_set = {len(resultant) for resultant in model.meta.exposure.read_pattern}

    @result_logger(requirement, 'Testing that the ramps are uneven')
    def test():
        assert len(length_set) > 1
    test()


@pytest.mark.bigdata
def test_science_verification(rampfit_result, rtdata_module, ignore_asdf_paths):
    """Check against expected data results"""
    requirement, _, _ = rampfit_result
    diff = compare_asdf(rtdata_module.output, rtdata_module.truth, **ignore_asdf_paths)

    @result_logger(requirement, 'Testing science veracity')
    def test():
        assert diff.identical, diff.report()
    test()


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
