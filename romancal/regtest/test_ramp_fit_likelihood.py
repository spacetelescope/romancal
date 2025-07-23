"""Regression test for the ramp fitting step with likelihood"""

import pytest
import roman_datamodels as rdm

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_liklihood_rampfit(
    rtdata, ignore_asdf_paths, resource_tracker, request, dms_logger
):
    """Testing ramp fitting  using the likelihood method"""

    input_data = "r0000101001001001001_0001_wfi01_f158_darkcurrent.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Define step (for running and log access)
    dms_logger.info("Testing ramp fitting with the likelihood algorithm. ")

    dms_logger.info(f"Image data file: {rtdata.input.rsplit('/', 1)[1]}")

    # Test Likelihood ramp fitting
    output = "r0000101001001001001_0001_wfi01_f158_like_rampfit.asdf"
    rtdata.output = output
    args = [
        "romancal.step.RampFitStep",
        rtdata.input,
        "--output_file=r0000101001001001001_0001_wfi01_f158_like_rampfit.asdf",
    ]
    dms_logger.info("Testing the likelihood fitting for ramps")
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    ramp_results = rdm.open(rtdata.output)

    dms_logger.info(
        "RampFit step recorded as complete? :"
        f" {ramp_results.meta.cal_step.ramp_fit == 'COMPLETE'}"
    )
    assert ramp_results.meta.cal_step.ramp_fit == "COMPLETE"

    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    dms_logger.info(f"Was the rate image data produced? : {diff.identical}")
    assert diff.identical, diff.report()
