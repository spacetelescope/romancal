from pathlib import Path

import pytest
import roman_datamodels as rdm

from romancal.lib.suffix import replace_suffix
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]

INPUT_FILENAME = "r0000101001001001001_0001_wfi01_f158_refpix.asdf"


@pytest.fixture(scope="module")
def run_dark_decay(rtdata_module, resource_tracker):
    rtdata = rtdata_module
    input_fn = INPUT_FILENAME
    output_fn = replace_suffix(Path(input_fn).stem, "dark_decay") + ".asdf"
    rtdata.get_data(f"WFI/image/{input_fn}")
    rtdata.output = output_fn
    rtdata.get_truth(f"truth/WFI/image/{output_fn}")
    with resource_tracker.track():
        RomanStep.from_cmdline(["romancal.step.DarkDecayStep", rtdata.input])
    return rtdata


@pytest.fixture(scope="module")
def output_model(run_dark_decay):
    with rdm.open(run_dark_decay.output) as model:
        yield model


def test_log_tracked_resources(log_tracked_resources, run_dark_decay):
    log_tracked_resources()


def test_output_matches_truth(run_dark_decay, ignore_asdf_paths):
    rtdata = run_dark_decay
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_cal_step_updated(output_model):
    """Test that cal_step is set as expected."""
    assert output_model.meta.cal_step.dark_decay == "COMPLETE"
