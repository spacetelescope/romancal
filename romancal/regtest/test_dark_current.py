"""Regression test for the Dark current subtraction step."""

import pytest

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_dark_current_subtraction_step(
    rtdata, ignore_asdf_paths, resource_tracker, request
):
    """Function to run and compare Dark Current subtraction files. Note: This
    should include tests for overrides etc."""

    input_datafile = "r0000101001001001001_0001_wfi01_linearity.asdf"
    rtdata.get_data(f"WFI/image/{input_datafile}")
    rtdata.input = input_datafile

    args = ["romancal.step.DarkCurrentStep", rtdata.input]
    with resource_tracker.track("dark"):
        RomanStep.from_cmdline(args)
    resource_tracker.log(request.node.properties)
    output = "r0000101001001001001_0001_wfi01_darkcurrent.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_dark_current_outfile_step(
    rtdata, ignore_asdf_paths, resource_tracker, request
):
    """Function to run and compare Dark Current subtraction files. Here the
    test is for renaming the output file."""
    input_datafile = "r0000101001001001001_0001_wfi01_linearity.asdf"
    rtdata.get_data(f"WFI/image/{input_datafile}")
    rtdata.input = input_datafile

    args = [
        "romancal.step.DarkCurrentStep",
        rtdata.input,
        "--output_file=Test_darkcurrent",
    ]
    with resource_tracker.track("dark"):
        RomanStep.from_cmdline(args)
    resource_tracker.log(request.node.properties)
    output = "Test_darkcurrent.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_dark_current_outfile_suffix(
    rtdata, ignore_asdf_paths, resource_tracker, request
):
    """Function to run and compare Dark Current subtraction files. Here the
    test is for renaming the output file."""
    input_datafile = "r0000101001001001001_0001_wfi01_linearity.asdf"
    rtdata.get_data(f"WFI/image/{input_datafile}")
    rtdata.input = input_datafile

    args = [
        "romancal.step.DarkCurrentStep",
        rtdata.input,
        "--output_file=Test_dark",
        '--suffix="suffix_test"',
    ]
    with resource_tracker.track("dark"):
        RomanStep.from_cmdline(args)
    resource_tracker.log(request.node.properties)
    output = "Test_darkcurrent.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_dark_current_output(rtdata, ignore_asdf_paths, resource_tracker, request):
    """Function to run and compare Dark Current subtraction files. Here the
    test for overriding the CRDS dark reference file."""

    input_datafile = "r0000101001001001001_0001_wfi01_linearity.asdf"
    rtdata.get_data(f"WFI/image/{input_datafile}")
    rtdata.input = input_datafile
    dark_output_name = "r0000101001001001001_0001_wfi01_darkcurrent.asdf"

    args = [
        "romancal.step.DarkCurrentStep",
        rtdata.input,
        f"--dark_output={dark_output_name}",
    ]
    with resource_tracker.track("dark"):
        RomanStep.from_cmdline(args)
    resource_tracker.log(request.node.properties)
    output = "r0000101001001001001_0001_wfi01_darkcurrent.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()
