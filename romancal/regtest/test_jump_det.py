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
    ]
    RomanStep.from_cmdline(args)
    output = "r0000101001001001001_01101_0001_WFI01_jump.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
@pytest.mark.parametrize(
    "base_name",
    [
        # "MLS0R1_alpha_1000_A_1000",
        "MLS0R1_alpha_1000_A_0",
        # "MLS0R1_alpha_1000_A_500",
        "MLS0R1_alpha_10_A_200",
        "MLS0R1_alpha_10_A_100",
        "MLS0R1_alpha_10_A_0",
    ],
)
def test_cas22_jump_detection(rtdata, ignore_asdf_paths, base_name):
    input_file = f"{base_name}_input.asdf"
    gain_file = f"{base_name}_gain.asdf"
    readnoise_file = f"{base_name}_readnoise.asdf"

    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.get_data(f"WFI/image/{gain_file}")
    rtdata.get_data(f"WFI/image/{readnoise_file}")
    rtdata.input = input_file
    rtdata.gain = gain_file
    rtdata.readnoise = readnoise_file

    args = [
        "romancal.step.RampFitStep",
        rtdata.input,
        f"--override_gain={rtdata.gain}",
        f"--override_readnoise={rtdata.readnoise}",
        "--use_jump_detection=True",
    ]
    RomanStep.from_cmdline(args)
    output = f"{base_name}_input_rampfit.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{base_name}_output.asdf")

    import asdf

    truth = asdf.open(rtdata.truth)["roman"]
    output = asdf.open(rtdata.output)["roman"]

    from numpy.testing import assert_allclose

    assert_allclose(truth["dq"], output["dq"])
    assert_allclose(truth["data"], output["data"])
    truth.close()
    output.close()
