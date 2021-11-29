import pytest

from romancal.stpipe import RomanStep
from .regtestdata import compare_asdf

@pytest.mark.bigdata
def test_ramp_fitting_step(rtdata, ignore_asdf_paths):
    rtdata.get_data("WFI/image/ramp.asdf")
    rtdata.input = "ramp.asdf"

    args = ["romancal.step.RampFitStep", rtdata.input, '--save_opt=True', '--opt_name=rampfit_opt.asdf']
    RomanStep.from_cmdline(args)
    output = "ramp_rampfitstep.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)

    output = "rampfit_opt_fitopt.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)
