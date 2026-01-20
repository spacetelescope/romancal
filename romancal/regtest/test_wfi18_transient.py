from pathlib import Path

import pytest
import roman_datamodels as rdm

from romancal.lib.suffix import replace_suffix
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]

WFI18_FILENAME = "TVAC2_NOMOPS_TTNOISE_20240418084515_WFI18_refpix.asdf"
WFI01_FILENAME = "r0000101001001001001_0001_wfi01_f158_refpix.asdf"


@pytest.fixture(
    params=[
        WFI18_FILENAME,
        WFI01_FILENAME,
    ],
    scope="module",
)
def run_wfi18_transient(rtdata_module, resource_tracker, request):
    rtdata = rtdata_module
    input_fn = request.param
    output_fn = replace_suffix(Path(input_fn).stem, "wfi18_transient") + ".asdf"
    rtdata.get_data(f"WFI/image/{input_fn}")
    rtdata.output = output_fn
    rtdata.get_truth(f"truth/WFI/image/{output_fn}")
    with resource_tracker.track():
        RomanStep.from_cmdline(["romancal.step.WFI18TransientStep", rtdata.input])
    return rtdata


@pytest.fixture(scope="module")
def output_model(run_wfi18_transient):
    with rdm.open(run_wfi18_transient.output) as model:
        yield model


def test_log_tracked_resources(log_tracked_resources, run_wfi18_transient):
    log_tracked_resources()


def test_output_matches_truth(run_wfi18_transient, ignore_asdf_paths):
    rtdata = run_wfi18_transient
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_cal_step_updated(output_model):
    """Test that cal_step is set as expected."""
    if output_model.meta.instrument.detector == "WFI18":
        assert output_model.meta.cal_step.wfi18_transient == "COMPLETE"
    else:
        assert output_model.meta.cal_step.wfi18_transient == "N/A"
