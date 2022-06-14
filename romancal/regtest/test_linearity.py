"""Regression test for the linearity correction step."""
import pytest

from romancal.stpipe import RomanStep
from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_linearity_step(rtdata, ignore_asdf_paths):
    """ Function to run and compare linearity correction files."""
    input_file = 'r0000101001001001001_01101_0001_WFI01_saturation.asdf'
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file

    args = ["romancal.step.LinearityStep", rtdata.input]
    RomanStep.from_cmdline(args)
    output =\
        "r0000101001001001001_01101_0001_WFI01_linearity.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth,
            **ignore_asdf_paths) is None)


@pytest.mark.bigdata
def test_linearity_outfile_step(rtdata, ignore_asdf_paths):
    """ Function to run and linearity correction files. Here the
        test is for renaming the output file. """
    input_file = 'r0000101001001001001_01101_0001_WFI01_saturation.asdf'
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file

    args = ["romancal.step.LinearityStep", rtdata.input,
            '--output_file=Test_linearity']
    RomanStep.from_cmdline(args)
    output = "Test_linearity.asdf"
    rtdata.output = output
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth,
            **ignore_asdf_paths) is None)
