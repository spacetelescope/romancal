"""Regression test for the linearity correction step."""

import pytest

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_linearity_step(rtdata, ignore_asdf_paths, resource_tracker, request):
    """Function to run and compare linearity correction files."""
    input_file = "r0000101001001001001_0001_wfi01_refpix.asdf"
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file

    args = ["romancal.step.LinearityStep", rtdata.input]
    with resource_tracker.track():
        RomanStep.from_cmdline(args)
    resource_tracker.log(request.node.user_properties)
    output = "r0000101001001001001_0001_wfi01_linearity.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_linearity_outfile_step(rtdata, ignore_asdf_paths, resource_tracker, request):
    """Function to run and linearity correction files. Here the
    test is for renaming the output file."""
    input_file = "r0000101001001001001_0001_wfi01_refpix.asdf"
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file

    args = ["romancal.step.LinearityStep", rtdata.input, "--output_file=Test_linearity"]
    with resource_tracker.track():
        RomanStep.from_cmdline(args)
    resource_tracker.log(request.node.user_properties)
    output = "Test_linearity.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()
