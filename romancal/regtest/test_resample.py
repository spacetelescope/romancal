import json
from io import StringIO

import numpy as np
import pytest
from metrics_logger.decorators import metrics_logger
from roman_datamodels import datamodels as rdm

from romancal.resample.resample_step import ResampleStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


def create_asn_file(
    output_filename: str = "resample_output.asdf",
    members_filename_list: list = None,
):
    asn_dict = {
        "asn_type": "None",
        "asn_rule": "DMS_ELPP_Base",
        "version_id": "null",
        "code_version": "0.9.1.dev28+ge987cc9.d20230106",
        "degraded_status": "No known degraded exposures in association.",
        "program": "noprogram",
        "constraints": "No constraints",
        "asn_id": "a3001",
        "target": "none",
        "asn_pool": "test_pool_name",
        "products": [
            {
                "name": output_filename,
                "members": [
                    {"expname": x, "exptype": "science"} for x in members_filename_list
                ],
            }
        ],
    }
    asn_content = json.dumps(asn_dict)
    asn_file_path = "sample_asn.json"
    asn_file = StringIO()
    asn_file.write(asn_content)
    with open(asn_file_path, mode="w") as f:
        print(asn_file.getvalue(), file=f)

    return asn_file_path


@metrics_logger(
    "DMS342"
)  # got DMS342 from here: https://jira.stsci.edu/browse/RSUBREQ-1051
@pytest.mark.bigdata
def test_resample_single_file(rtdata, ignore_asdf_paths):
    input_data = [
        "r0000501001001001001_01101_0001_WFI02_cal_proc_resample.asdf",
        "r0000501001001001001_01101_0002_WFI02_cal_proc_resample.asdf",
    ]
    output_data = "resample_output_resamplestep.asdf"
    truth_data = "resample_truth_resamplestep.asdf"

    [rtdata.get_data(f"WFI/image/{data}") for data in input_data]
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = create_asn_file(members_filename_list=input_data)
    rtdata.output = output_data

    # instantiate ResampleStep (for running and log access)
    step = ResampleStep()

    args = [
        "romancal.step.ResampleStep",
        rtdata.input,
        "--rotation=0",
        f"--output_file='{rtdata.output}'",
    ]
    RomanStep.from_cmdline(args)
    resample_out = rdm.open(rtdata.output)

    step.log.info(
        "ResampleStep recorded as complete? :"
        f' {resample_out.meta.cal_step.resample == "COMPLETE"}'
    )
    assert resample_out.meta.cal_step.resample == "COMPLETE"

    step.log.info(
        "ResampleStep created 'meta.resample'? :"
        f' {hasattr(resample_out.meta, "resample")}'
    )
    assert hasattr(resample_out.meta, "resample")

    step.log.info(
        f"""DMS343 MSG: ResampleStep created new attribute data quality information? :\
            {
                all(
                    hasattr(resample_out, x) for x in [
                        "data",
                        "err",
                        "var_poisson",
                        "var_rnoise",
                        "var_flat",
                    ]
                )
            }"""
    )
    assert all(
        hasattr(resample_out, x)
        for x in ["data", "err", "var_poisson", "var_rnoise", "var_flat"]
    )

    step.log.info(
        f"""DMS343 MSG: Were the variance arrays populated (variance propagation)? :\
            {
                all(
                    np.sum(~np.isnan(getattr(resample_out, x))) for x in [
                        "var_poisson",
                        "var_rnoise",
                    ]
                )
            }"""  # noqa: E501
    )
    assert all(
        np.sum(~np.isnan(getattr(resample_out, x)))
        for x in ["var_poisson", "var_rnoise"]
    )

    step.log.info(
        f"""DMS343 MSG: Are there NaNs or zeros in the variance arrays, indicating poor data quality? :\
            {
                any(
                    np.sum(
                        np.logical_or(
                            np.isnan(getattr(resample_out, x)),
                            np.equal(getattr(resample_out, x), 0)
                        )
                    ) > 0 for x in ["var_poisson", "var_rnoise", "var_flat"]
                )

            }"""
    )
    assert all(
        np.sum(np.isnan(getattr(resample_out, x)))
        for x in ["var_poisson", "var_rnoise", "var_flat"]
    )

    step.log.info(
        f"""DMS344 MSG: ResampleStep created new attribute with total exposure time? :\
            {"product_exposure_time" in resample_out.meta.resample}
        """
    )
    assert "product_exposure_time" in resample_out.meta.resample

    step.log.info(
        f"""DMS345 MSG: ResampleStep included all metadata relevant to the creation of the mosaic? :\
            {
                all(
                    hasattr(resample_out.meta.resample, x)
                    and bool(getattr(resample_out.meta.resample, x))
                    for x in [
                        "pixel_scale_ratio",
                        "pixfrac",
                        "pointings",
                        "product_exposure_time",
                        "weight_type",
                        "members",
                    ]
                )
            }"""  # noqa: E501
    )
    assert all(
        hasattr(resample_out.meta.resample, x)
        and bool(getattr(resample_out.meta.resample, x))
        for x in [
            "pixel_scale_ratio",
            "pixfrac",
            "pointings",
            "product_exposure_time",
            "weight_type",
            "members",
        ]
    )

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    step.log.info("Was the proper Resample data produced?" f" : {diff.identical}")
    assert diff.identical, diff.report()
