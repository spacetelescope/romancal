"""Regression test for the Jump detection step."""
import pytest

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_jump_detection_step(rtdata, ignore_asdf_paths):
    """Function to run and compare Jump Detection files. Note: This should
    include tests for overrides etc."""

    input_file = "r0000101001001001001_01101_0001_WFI01_darkcurrent.asdf"
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file

    # Note: the thresholds should be reset to the defaults once we have better
    # input data
    args = [
        "romancal.step.JumpStep",
        rtdata.input,
        "--rejection_threshold=180.",
        "--three_group_rejection_threshold=185.",
        "--four_group_rejection_threshold=190.",
        # disable the jump detection in the rampfit step so this one will run
        "--run_ramp_jump_detection=False",
    ]
    RomanStep.from_cmdline(args)
    output = "r0000101001001001001_01101_0001_WFI01_jump.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_cas22_jump_detection(rtdata, ignore_asdf_paths):
    input_file = "r0000101001001001001_01101_0001_WFI01_darkcurrent.asdf"
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file
    args = [
        "romancal.step.RampFitStep",
        rtdata.input,
        # Ensure the jump detection is run in the rampfit step
        "--use_ramp_jump_detection=True",
        # Add differentiating suffix to output file (so it doesn't collide with rampfit)
        "--suffix=rampfit_jump",
    ]
    RomanStep.from_cmdline(args)
    output = "r0000101001001001001_01101_0001_WFI01_rampfit_jump.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()
