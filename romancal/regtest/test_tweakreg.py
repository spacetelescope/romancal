import pytest
from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep
from romancal.tweakreg.tweakreg_step import TweakRegStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_is_wcs_correction_small(rtdata, ignore_asdf_paths):
    input_data = "r0000401001001001001_01101_0001_WFI01_cal_tweakreg_input.asdf"
    output_data = "r0000401001001001001_01101_0001_WFI01_cal_tweakreg_output.asdf"
    truth_data = "r0000401001001001001_01101_0001_WFI01_cal_tweakreg_truth.asdf"

    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = input_data
    rtdata.output = output_data

    step = TweakRegStep()

    args = [
        "romancal.step.TweakRegStep",
        [rtdata.input],
        f"--output_file='{rtdata.output}'",
    ]
    RomanStep.from_cmdline(args)
    tweakreg_out = rdm.open(rtdata.output)

    step.log.info(
        "DMSXXX MSG: TweakReg step recorded as complete? :"
        f' {tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"}'
    )
    assert tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"

    step.log.info(
        "DMSXXX MSG: Was the proper TweakReg data produced?"
        f" : {(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}"
    )
    assert compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None
