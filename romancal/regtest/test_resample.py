import numpy as np
import pytest
from metrics_logger.decorators import metrics_logger
from roman_datamodels import datamodels as rdm

from romancal.resample.resample_step import ResampleStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@metrics_logger("DMS342", "DMS343", "DMS344", "DMS345")
@pytest.mark.bigdata
def test_resample_single_file(rtdata, ignore_asdf_paths):
    input_data = [
        "r0000101001001001001_01101_0001_WFI01_cal.asdf",
        "r0000101001001001001_01101_0002_WFI01_cal.asdf",
    ]
    output_data = "mosaic_resamplestep.asdf"

    [rtdata.get_data(f"WFI/image/{data}") for data in input_data]
    asnfn = "mosaic_asn.json"
    rtdata.get_data(f"WFI/image/{asnfn}")
    rtdata.get_truth(f"truth/WFI/image/{output_data}")

    rtdata.input = asnfn
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
        f"""DMS342 MSG: Was ICRS used as the mosaic astrometric reference frame? :\
            {
                resample_out.meta.coordinates.reference_frame == "ICRS"
            }
        """
    )
    assert resample_out.meta.coordinates.reference_frame == "ICRS"

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
