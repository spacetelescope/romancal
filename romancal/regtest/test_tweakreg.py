import pytest
from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep
from romancal.tweakreg.tweakreg_step import TweakRegStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_tweakreg(rtdata, ignore_asdf_paths):
    # N.B.: the data were created using WFIsim and processed through
    # the three pipeline steps listed below:
    # - assign_wcs;
    # - photom;
    # - source_detection.
    input_data = "r0000401001001001001_01101_0001_WFI01_cal_tweakreg.asdf"
    output_data = "r0000401001001001001_01101_0001_WFI01_output.asdf"
    truth_data = "r0000401001001001001_01101_0001_WFI01_cal_twkreg_proc.asdf"

    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = input_data
    rtdata.output = output_data

    # add filename (this is set in exposure_pipeline right after running dq_init)
    # with rdm.open(rtdata.input) as model:
    #     model.meta["filename"] = Path(rtdata.input).name
    #     model.save(Path(rtdata.input))

    # instantiate TweakRegStep (for running and log access)
    step = TweakRegStep()

    args = [
        "romancal.step.TweakRegStep",
        [rtdata.input],
        f"--output_file='{rtdata.output}'",
        "--suffix='output'",
    ]
    RomanStep.from_cmdline(args)
    tweakreg_out = rdm.open(rtdata.output)

    step.log.info(
        "DMS280 MSG: TweakReg step recorded as complete? :"
        f' {tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"}'
    )
    assert tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"

    step.log.info(
        f"""DMS280 MSG: TweakReg created new attribute with fit results? :\
            {"wcs_fit_results" in tweakreg_out.meta}"""
    )
    assert "wcs_fit_results" in tweakreg_out.meta

    step.log.info(
        f"""DMS280 MSG: TweakReg created new coordinate frame 'v2v3corr'? :\
            {"v2v3corr" in tweakreg_out.meta.wcs.available_frames}"""
    )
    assert "v2v3corr" in tweakreg_out.meta.wcs.available_frames

    step.log.info(
        "DMS280 MSG: Was the proper TweakReg data produced?"
        f" : {(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}"
    )
    assert compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None
