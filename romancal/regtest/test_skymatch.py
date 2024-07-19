import pytest
from romancal.stpipe import RomanStep
from roman_datamodels import datamodels as rdm

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_skymatch_step(rtdata, ignore_asdf_paths):
    """Test for the flat field step using imaging data."""

    l2_data_folder = "/Users/mteodoro/ROMAN/SYNTHETIC_IMAGES/IMAGES/EDDIE/romanisim/20240311/PROCESSED/ELP"

    input_data = [
        "r0000101001001001001_01101_0001_WFI01_cal.asdf",
        "r0000101001001001001_01101_0002_WFI01_cal.asdf",
    ]
    output_data = "skymatchstep.asdf"

    [rtdata.get_data(f"{l2_data_folder}/{data}") for data in input_data]
    asnfn = "skymatch_asn.json"
    rtdata.get_data(f"{l2_data_folder}/{asnfn}")
    # rtdata.get_truth(f"truth/WFI/image/{output_data}")

    rtdata.input = asnfn
    rtdata.output = output_data

    args = [
        "romancal.step.SkyMatchStep",
        rtdata.input,
        f"--output_file='{rtdata.output}'",
    ]
    RomanStep.from_cmdline(args)

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)

    assert diff.identical, diff.report()
