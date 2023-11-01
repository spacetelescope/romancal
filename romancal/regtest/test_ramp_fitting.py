""" Module to test rampfit with optional output """
import functools
from pathlib import Path
import pytest

import roman_datamodels as rdm
from romancal.lib.suffix import replace_suffix
from romancal.ramp_fitting import RampFitStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf

# Logger to use for test result reporting
logger = RampFitStep().log


def result_logger(message, logger=logger):
    """Decorator to log assertion status

    Parameters
    ----------
    message : str
        The message to use. Note that the phrase of either "PASS" or "FAIL"
        will be appended to the end of the message.
    """
    def decorator(test_function):

        @functools.wraps(test_function)
        def inner(*args, **kwargs):

            def log_result(result):
                result_text = 'PASS' if result else 'FAIL'
                logger.info(message + result_text)

            try:
                test_function(*args, **kwargs)
            except Exception:
                log_result(False)
                raise
            else:
                log_result(True)

        return inner

    return decorator


@pytest.mark.bigdata
@result_logger('DMS362 MSG: Testing that the result file is of type "asdf".......')
def test_is_asdf(rampfit_result):
    """Check that the filename has the correct file type"""
    _, result_path = rampfit_result

    assert result_path.exists() and result_path.suffix == '.asdf'


@pytest.mark.bigdata
@result_logger('DMS362 MSG: Testing that the result model is Level 2.......')
def test_is_imagemodel(rampfit_result):
    """Check that the result is an ImageModel"""
    model, _ = rampfit_result

    assert isinstance(model, rdm.datamodels.ImageModel)


@pytest.mark.bigdata
@result_logger('DMS362 MSG: Testing that the result file has the suffix "rampfit".......')
def test_is_rampfit(rampfit_result):
    """Check that the calibration suffix is 'rampfit'"""
    _, result_path = rampfit_result

    assert result_path.exists() and result_path.stem.endswith('rampfit')


@pytest.mark.bigdata
@result_logger('DMS362 MSG: Testing that RampFitStep completed.......')
def test_is_step_complete(rampfit_result):
    """Check that the calibration step is marked complete"""
    model, _ = rampfit_result

    assert model.meta.cal_step.ramp_fit == 'COMPLETE'


@pytest.mark.bigdata
@result_logger('DMS362 MSG: Testing that the ramps are uneven.......')
def test_is_uneven(rampfit_result):
    """Verify that the provided model represents uneven ramps

    Parameters
    ----------
    rampfit_result : `roman_datamodels.ImageModel`
        Model created from `RampFitStep`
    """
    model, _ = rampfit_result
    length_set = {len(resultant) for resultant in model.meta.exposure.read_pattern}

    assert len(length_set) > 1


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


# ######################
# fixtures and utilities
# ######################
@pytest.fixture(scope='module')
def rampfit_result(rtdata_module):
    """Run RampFitStep

    Parameters
    ----------
    rtdata_module : pytest.fixture
        artifactory fixture for data retrieval

    Returns
    -------
    model, path, logger : `ImageModel`, `pathlib.Path`
        Model and path to model.
        Logger is used to log statuses of the tests being run.
    """
    input_data = 'random_dqinit.asdf'
    input_data = rtdata_module.get_data(f'WFI/image/{input_data}')
    result_model = RampFitStep.call(input_data, save_results=True)

    input_data_path = Path(input_data)
    expected_path = input_data_path.parent / (replace_suffix(input_data_path.stem, 'rampfit') + '.asdf')

    try:
        yield result_model, expected_path
    finally:
        result_model.close()


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"
