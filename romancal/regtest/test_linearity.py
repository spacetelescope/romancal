"""Regression test for the linearity correction step."""
import pytest

from romancal.stpipe import RomanStep
from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_linearity_step(rtdata, ignore_asdf_paths):
    """ Function to run and compare linearity correction files."""
    rtdata.get_data(
        "WFI/image/r0000101001001001001_01101_0001_WFI01_saturationstep.asdf")
    rtdata.input = "r0000101001001001001_01101_0001_WFI01_saturationstep.asdf"

    args = ["romancal.step.LinearityStep", rtdata.input]
    RomanStep.from_cmdline(args)
    output =\
        "r0000101001001001001_01101_0001_WFI01_linearitystep.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth,
            **ignore_asdf_paths) is None)


@pytest.mark.bigdata
def test_linearity_outfile_step(rtdata, ignore_asdf_paths):
    """ Function to run and linearity correction files. Here the
        test is for renaming the output file. """
    rtdata.get_data(
        "WFI/image/r0000101001001001001_01101_0001_WFI01_saturationstep.asdf")
    rtdata.input = "r0000101001001001001_01101_0001_WFI01_saturationstep.asdf"

    args = ["romancal.step.LinearityStep", rtdata.input,
            '--output_file=Test_linearitystep']
    RomanStep.from_cmdline(args)
    output = "Test_linearitystep.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth,
            **ignore_asdf_paths) is None)
