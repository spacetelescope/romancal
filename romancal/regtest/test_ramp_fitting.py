""" Module to test rampfit with optional output """
from pathlib import Path
import pytest

import roman_datamodels as rdm
from romancal.lib.suffix import replace_suffix
from romancal.ramp_fitting import RampFitStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_is_imagemodel(rampfit_result):
    """Check that the result is an ImageModel"""
    model, _ = rampfit_result

    assert isinstance(model, rdm.datamodels.ImageModel)


@pytest.mark.bigdata
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
    model, path : `ImageModel`, `pathlib.Path`
        Model and path to model
    """
    input_data = 'random_dqinit.asdf'
    input_data = rtdata_module.get_data(f'WFI/image/{input_data}')
    with rdm.open(input_data) as input_model:
        result_model = RampFitStep.call(input_model, save_results=True)

    expected_path = replace_suffix(Path(input_data).stem, 'rampfit') + '.asdf'

    try:
        yield result_model, expected_path
    finally:
        result_model.close()
