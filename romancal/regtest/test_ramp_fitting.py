""" Module to test rampfit with optional output """

import pytest
from romancal.stpipe import RomanStep
from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_ramp_fitting_step(rtdata, ignore_asdf_paths):
    """ Testing the ramp fitting step"""
    input_data = "r0000101001001001001_01101_0001_WFI01_jump.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    args = ["romancal.step.RampFitStep", rtdata.input, '--save_opt=True',
            '--opt_name=rampfit_opt.asdf']
    RomanStep.from_cmdline(args)
    output = "r0000101001001001001_01101_0001_WFI01_rampfit.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None

    output = "rampfit_opt_fitopt.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert compare_asdf(rtdata.output, rtdata.truth,
                        **ignore_asdf_paths) is None
