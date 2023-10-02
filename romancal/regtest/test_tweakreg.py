from io import StringIO

import pytest
from metrics_logger.decorators import metrics_logger
from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep
from romancal.tweakreg.tweakreg_step import TweakRegStep

from .regtestdata import compare_asdf


def create_asn_file():
    asn_content = """
        {
            "asn_type": "None",
            "asn_rule": "DMS_ELPP_Base",
            "version_id": null,
            "code_version": "0.9.1.dev28+ge987cc9.d20230106",
            "degraded_status": "No known degraded exposures in association.",
            "program": "noprogram",
            "constraints": "No constraints",
            "asn_id": "a3001",
            "target": "none",
            "asn_pool": "test_pool_name",
            "products": [
                {
                    "name": "files.asdf",
                    "members": [
                        {
                            "expname": "r0000501001001001001_01101_0001_WFI02_cal_tweakreg.asdf",
                            "exptype": "science"
                        }
                    ]
                }
            ]
        }
    """  # noqa: E501
    asn_file_path = "sample_asn.json"
    asn_file = StringIO()
    asn_file.write(asn_content)
    with open(asn_file_path, mode="w") as f:
        print(asn_file.getvalue(), file=f)

    return asn_file_path


@metrics_logger("DMS280")
@pytest.mark.bigdata
def test_tweakreg(rtdata, ignore_asdf_paths, tmp_path):
    # N.B.: the data were created using WFIsim and processed through
    # the three pipeline steps listed below:
    # - assign_wcs;
    # - photom;
    # - source_detection.
    input_data = "r0000501001001001001_01101_0001_WFI02_cal_tweakreg.asdf"
    output_data = "r0000501001001001001_01101_0001_WFI02_output.asdf"
    truth_data = "r0000501001001001001_01101_0001_WFI02_tweakreg.asdf"

    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = create_asn_file()
    rtdata.output = output_data

    # instantiate TweakRegStep (for running and log access)
    step = TweakRegStep()

    args = [
        "romancal.step.TweakRegStep",
        rtdata.input,
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

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    step.log.info(
        f"DMS280 MSG: Was the proper TweakReg data produced? : {diff.identical}"
    )
    assert diff.identical, diff.report()
