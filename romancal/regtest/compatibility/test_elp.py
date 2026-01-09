from pathlib import Path

import pytest

from romancal.pipeline.exposure_pipeline import ExposurePipeline

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_elp(old_rtdata_module, resource_tracker):
    rtdata = old_rtdata_module

    input_data = "r0000101001001001001_0001_wfi01_f158_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000101001001001001_0001_wfi01_f158_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
        "--update_version=True",
    ]
    with resource_tracker.track():
        ExposurePipeline.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_elp):
    return run_elp.output


def test_log_tracked_resources(log_tracked_resources, run_elp):
    log_tracked_resources()


def test_output_exists(output_filename):
    assert Path(output_filename).exists()
