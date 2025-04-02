"""Roman tests for flat field correction"""

import pytest
import roman_datamodels as rdm

from romancal.pipeline.exposure_pipeline import ExposurePipeline

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_elp(rtdata_module, resource_tracker):
    rtdata = rtdata_module

    # The input data is from INS for stress testing at some point this should be generated
    # by INS every time new data is needed.

    input_data = "r0000201001001001001_0004_wfi01_grism_uncal.asdf"
    rtdata.get_data(f"WFI/grism/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000201001001001001_0004_wfi01_grism_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
    ]
    with resource_tracker.track():
        ExposurePipeline.from_cmdline(args)
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_elp):
    return run_elp.output


@pytest.fixture(scope="module")
def output_model(output_filename):
    with rdm.open(output_filename) as model:
        yield model


@pytest.fixture(scope="module")
def input_filename(run_elp):
    return run_elp.input


@pytest.fixture(scope="module")
def input_model(input_filename):
    with rdm.open(input_filename) as model:
        yield model


def test_log_tracked_resources(log_tracked_resources, run_elp):
    log_tracked_resources()


def test_output_is_image_model(output_model):
    # DMS414
    assert isinstance(output_model, rdm.datamodels.ImageModel)


def test_input_has_16_resultants(input_model):
    # DMS414
    assert len(input_model.meta.exposure.read_pattern) == 16


def test_output_has_16_resultants(output_model):
    # DMS414
    assert len(output_model.meta.exposure.read_pattern) == 16


@pytest.mark.parametrize(
    "step_name, status",
    (
        ("assign_wcs", "COMPLETE"),
        ("dark", "COMPLETE"),
        ("dq_init", "COMPLETE"),
        ("linearity", "COMPLETE"),
        ("ramp_fit", "COMPLETE"),
        ("saturation", "COMPLETE"),
        # skipped
        ("flat_field", "SKIPPED"),
        ("photom", "SKIPPED"),
    ),
)
def test_step_status(output_model, step_name, status):
    assert getattr(output_model.meta.cal_step, step_name) == status
