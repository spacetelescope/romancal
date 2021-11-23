"""Regression test for the Jump detection step."""
import pytest

from romancal.stpipe import RomanStep
from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_jump_detection_step(rtdata, ignore_asdf_paths):
    """ Function to run and compare Jump Detection files. Note: This should
        include tests for overrides etc. """
    rtdata.get_data("WFI/image/l1_0007_science_raw_dqinitstep.asdf")
    rtdata.input = "l1_0007_science_raw_dqinitstep.asdf"

    args = ["romancal.step.JumpStep", rtdata.input]
    RomanStep.from_cmdline(args)
    output = "l1_0007_science_raw_dqinitstep_jumpstep.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth,
            **ignore_asdf_paths) is None)
