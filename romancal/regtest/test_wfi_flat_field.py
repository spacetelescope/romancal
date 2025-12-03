from pathlib import Path

import pytest
import roman_datamodels as rdm
from astropy.time import Time

from romancal.lib.suffix import replace_suffix
from romancal.step import FlatFieldStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]

IMAGE_FILENAME = "r0000101001001001001_0001_wfi01_f158_assignwcs.asdf"
GRISM_FILENAME = "r0000201001001001001_0001_wfi01_grism_assignwcs.asdf"


@pytest.fixture(
    params=[
        IMAGE_FILENAME,
        GRISM_FILENAME,
    ],
    scope="module",
)
def run_flat_field(rtdata_module, resource_tracker, request):
    rtdata = rtdata_module
    input_fn = request.param
    output_fn = replace_suffix(Path(input_fn).stem, "flat") + ".asdf"
    file_type = "grism" if "grism" in input_fn else "image"
    rtdata.get_data(f"WFI/{file_type}/{input_fn}")
    rtdata.output = output_fn
    rtdata.get_truth(f"truth/WFI/{file_type}/{output_fn}")
    with resource_tracker.track():
        RomanStep.from_cmdline(["romancal.step.FlatFieldStep", rtdata.input])
    return rtdata


@pytest.fixture(scope="module")
def output_model(run_flat_field):
    with rdm.open(run_flat_field.output) as model:
        yield model


def test_log_tracked_resources(log_tracked_resources, run_flat_field):
    log_tracked_resources()


def test_output_matches_truth(run_flat_field, ignore_asdf_paths):
    rtdata = run_flat_field
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_cal_step_updated(output_model):
    """Test that cal_step is step as expected"""
    if output_model.meta.exposure.type == "WFI_IMAGE":
        assert output_model.meta.cal_step.flat_field == "COMPLETE"
    else:
        assert output_model.meta.cal_step.flat_field == "SKIPPED"


def test_ref_file_used(output_model):
    """Test that a ref file was used and recorded where expected"""
    if output_model.meta.exposure.type == "WFI_IMAGE":
        assert "roman_wfi_flat" in output_model.meta.ref_file.flat
    else:
        assert output_model.meta.ref_file.flat == "N/A"


@pytest.mark.xfail(
    reason="Cannot succeed because no references yet exist with more than one USEAFTER"
)
def test_ref_file_query(rtdata):
    """Test DMS79: 2 different start times return different flats"""
    rtdata.get_data(f"WFI/image/{IMAGE_FILENAME}")
    step = FlatFieldStep()
    with rdm.open(IMAGE_FILENAME) as model:
        ref_a = step.get_reference_file(model, "flat")
        model.meta.exposure.start_time = Time("2020-01-01T00:00:00", format="isot")
        ref_b = step.get_reference_file(model, "flat")
    assert ref_a != "N/A"
    assert ref_b != "N/A"
    assert ref_a != ref_b
