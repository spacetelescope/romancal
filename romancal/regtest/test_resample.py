import numpy as np
import pytest
from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_resample_single_file(
    rtdata, ignore_asdf_paths, resource_tracker, request, dms_logger
):
    output_data = "mosaic_resamplestep.asdf"

    rtdata.get_asn("WFI/image/L3_mosaic_asn.json")
    rtdata.get_truth(f"truth/WFI/image/{output_data}")

    rtdata.output = output_data

    args = [
        "romancal.step.ResampleStep",
        rtdata.input,
        "--rotation=0",
        "--resample_on_skycell=False",
        f"--output_file='{rtdata.output}'",
    ]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    resample_out = rdm.open(rtdata.output)

    dms_logger.info(
        "ResampleStep recorded as complete? :"
        f" {resample_out.meta.cal_step.resample == 'COMPLETE'}"
    )
    assert resample_out.meta.cal_step.resample == "COMPLETE"

    dms_logger.info(
        "ResampleStep created 'meta.resample'? :"
        f" {hasattr(resample_out.meta, 'resample')}"
    )
    assert hasattr(resample_out.meta, "resample")

    dms_logger.info(
        f"""DMS342 MSG: Was ICRS used as the mosaic astrometric reference frame? :\
            {resample_out.meta.coordinates.reference_frame == "ICRS"}
        """
    )
    assert resample_out.meta.coordinates.reference_frame == "ICRS"

    dms_logger.info(
        f"""DMS343 MSG: ResampleStep created new attribute data quality information? :\
            {
            all(
                hasattr(resample_out, x)
                for x in [
                    "data",
                    "err",
                    "var_poisson",
                    "var_rnoise",
                ]
            )
        }"""
    )
    assert all(
        hasattr(resample_out, x) for x in ["data", "err", "var_poisson", "var_rnoise"]
    )

    dms_logger.info(
        f"""DMS343 MSG: Were the variance arrays populated (variance propagation)? :\
            {
            all(
                np.sum(~np.isnan(getattr(resample_out, x)))
                for x in [
                    "var_poisson",
                    "var_rnoise",
                ]
            )
        }"""
    )
    assert all(
        np.sum(~np.isnan(getattr(resample_out, x)))
        for x in ["var_poisson", "var_rnoise"]
    )

    dms_logger.info(
        f"""DMS343 MSG: Are there NaNs or zeros in the variance arrays, indicating poor data quality? :\
            {
            any(
                np.sum(
                    np.logical_or(
                        np.isnan(getattr(resample_out, x)),
                        np.equal(getattr(resample_out, x), 0),
                    )
                )
                > 0
                for x in ["var_poisson", "var_rnoise"]
            )
        }"""
    )
    assert all(
        np.sum(np.isnan(getattr(resample_out, x)))
        for x in ["var_poisson", "var_rnoise"]
    )

    dms_logger.info(
        f"""DMS344 MSG: ResampleStep created new attribute with total exposure time? :\
            {"max_exposure_time" in resample_out.meta.coadd_info}
        """
    )
    assert "max_exposure_time" in resample_out.meta.coadd_info

    dms_logger.info(
        f"""DMS345 MSG: ResampleStep included all metadata relevant to the creation of the mosaic? :\
            {
            all(
                hasattr(resample_out.meta.resample, x)
                and bool(getattr(resample_out.meta.resample, x))
                for x in [
                    "pixel_scale_ratio",
                    "pixfrac",
                    "pointings",
                    "weight_type",
                    "members",
                ]
            )
        }"""
    )
    assert all(
        hasattr(resample_out.meta.resample, x)
        and bool(getattr(resample_out.meta.resample, x))
        for x in [
            "pixel_scale_ratio",
            "pixfrac",
            "pointings",
            "weight_type",
            "members",
        ]
    )

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    dms_logger.info(f"Was the proper Resample data produced? : {diff.identical}")
    assert diff.identical, diff.report()
