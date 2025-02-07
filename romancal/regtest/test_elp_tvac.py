"""Test that TVAC/FPS run in the pipeline"""

import numpy as np
import pytest
import roman_datamodels as rdm
from gwcs.wcstools import grid_from_bounding_box
from roman_datamodels.dqflags import pixel

from romancal.pipeline.exposure_pipeline import ExposurePipeline

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = pytest.mark.bigdata


@pytest.fixture(scope="module")
def run_elp(rtdata_module):
    rtdata = rtdata_module

    # Get reference
    rtdata.get_data('references/dark_ma510.asdf')

    input_data = "TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
        "--steps.dark_current.override_dark=dark_ma510.asdf",
        "--steps.rampfit.override_dark=dark_ma510.asdf"
    ]
    ExposurePipeline.from_cmdline(args)

    # get truth file
    rtdata.get_truth(f"truth/WFI/image/{output}")
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_elp):
    return run_elp.output


@pytest.fixture(scope="module")
def output_model(output_filename):
    with rdm.open(output_filename) as model:
        yield model


@pytest.fixture(scope="module")
def truth_filename(run_elp):
    return run_elp.truth


@pytest.mark.soctests
def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.soctests
def test_output_is_image_model(output_model):
    # DMS280 result is an ImageModel
    assert isinstance(output_model, rdm.datamodels.ImageModel)


@pytest.mark.soctests
@pytest.mark.parametrize(
    "step_name",
    (
        "assign_wcs",
        "flat_field",
        "dark",
        "dq_init",
        "linearity",
        "ramp_fit",
        "saturation",
    ),
)
def test_steps_ran(output_model, step_name):
    # DMS86
    # also DMS129 for assign_wcs
    assert getattr(output_model.meta.cal_step, step_name) == "COMPLETE"


@pytest.mark.soctests
def test_has_a_wcs(output_model):
    # DMS129
    assert output_model.meta.wcs is not None


@pytest.mark.soctests
def test_wcs_has_distortion_information(output_model):
    # DMS129
    assert "v2v3" in output_model.meta.wcs.available_frames


@pytest.mark.soctests
def test_wcs_applies_distortion_correction(output_model):
    # DMS129
    # compare coordinates before and after distortion correction has been applied
    # 1 - get new image array based on the model
    x0, y0 = grid_from_bounding_box(output_model.meta.wcs.bounding_box)
    # 2 - apply the distortion-corrected WCS solution to new image array
    corrected_coords = output_model.meta.wcs(x0, y0)
    # 3 - apply the transformation from 'v2v3' to 'world' without distortion correction
    original_coords = output_model.meta.wcs.get_transform("v2v3", "world")(x0, y0)
    # compare both results to make sure they don't match
    # (which means the distortion correction was actually applied to the model)
    assert (corrected_coords[0] != original_coords[0]).all()
    assert (corrected_coords[1] != original_coords[1]).all()


@pytest.mark.soctests
@pytest.mark.parametrize(
    "arr_name", ("dq", "err", "var_poisson", "var_rnoise", "var_flat")
)
def test_array_exists(output_model, arr_name):
    # DMS87
    assert hasattr(output_model, arr_name)


@pytest.mark.soctests
def test_has_exposure_time(output_model):
    # DMS88 total exposure time exists
    assert "exposure_time" in output_model.meta.exposure


@pytest.mark.soctests
@pytest.mark.parametrize("meta_attribute", ("detector", "optical_element"))
def test_instrument_meta(output_model, meta_attribute):
    # DMS136 PSF tests
    assert meta_attribute in output_model.meta.instrument


@pytest.mark.soctests
def test_wcs_has_bounding_box(output_model):
    # DMS89 WCS tests
    assert len(output_model.meta.wcs.bounding_box) == 2


def test_elp_input_dm(rtdata, ignore_asdf_paths):
    """Test for input roman Datamodel to exposure level pipeline"""
    input_data = "r0000101001001001001_0001_wfi01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    dm_input = rdm.open(rtdata.input)

    # Test Pipeline with input datamodel
    output = "r0000101001001001001_0001_wfi01_cal.asdf"
    rtdata.output = output
    ExposurePipeline.call(dm_input, save_results=True)
    rtdata.get_truth(f"truth/WFI/image/{output}")

    # Ensure step completion is as expected
    with rdm.open(rtdata.output) as model:
        assert model.meta.cal_step.dq_init == "COMPLETE"
        assert model.meta.cal_step.saturation == "COMPLETE"
        assert model.meta.cal_step.linearity == "COMPLETE"
        assert model.meta.cal_step.dark == "COMPLETE"
        assert model.meta.cal_step.ramp_fit == "COMPLETE"
        assert model.meta.cal_step.assign_wcs == "COMPLETE"
        assert model.meta.cal_step.flat_field == "COMPLETE"
        assert model.meta.cal_step.photom == "COMPLETE"
