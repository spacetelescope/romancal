import pytest

from romancal.pipeline.exposure_pipeline import ExposurePipeline

from ..regtestdata import compare_asdf

# mark all tests in this module
pytestmark = pytest.mark.bigdata


@pytest.fixture(scope="module")
def run_elp(rtdata_module, resource_tracker):
    rtdata = rtdata_module

    input_data = "r0000101001001001001_0001_wfi01_f158_uncal.asdf"
    rtdata.get_data(f"WFI/image/compatibility/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000101001001001001_0001_wfi01_f158_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
    ]
    with resource_tracker.track():
        ExposurePipeline.from_cmdline(args)

    # get truth file
    rtdata.get_truth(f"truth/WFI/image/compatbility/{output}")
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_elp):
    return run_elp.output


@pytest.fixture(scope="module")
def truth_filename(run_elp):
    return run_elp.truth


def test_log_tracked_resources(log_tracked_resources, run_elp):
    log_tracked_resources()


@pytest.mark.soctests
def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()
