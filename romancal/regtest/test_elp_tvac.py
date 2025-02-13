"""Test that TVAC/FPS run in the pipeline"""

from pathlib import Path
import pytest
import roman_datamodels as rdm
from gwcs.wcstools import grid_from_bounding_box

from romancal.lib.suffix import replace_suffix
from romancal.pipeline.exposure_pipeline import ExposurePipeline

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = pytest.mark.bigdata


@pytest.fixture(scope="module")
def run_elp(rtdata_module):
    rtdata = rtdata_module

    # Get reference
    rtdata.get_data("references/dark_ma510.asdf")

    input_data = "TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_cal.asdf"
    output_dqinit = "TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_dqinit.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
        "--steps.dq_init.save=true",
        "--steps.dark_current.override_dark=dark_ma510.asdf",
        "--steps.rampfit.override_dark=dark_ma510.asdf",
    ]
    ExposurePipeline.from_cmdline(args)

    # get truth file
    rtdata.get_truth(f"truth/WFI/image/{output_dqinit}")
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


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_dqinit_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    dqinit_path = Path(output_filename)
    dqinit_filename = replace_suffix(dqinit_path.stem, 'dqinit')
    dqinit_filename = dqinit_path.parent / (dqinit_filename + dqinit_path.suffix)

    truth_path = Path(truth_filename)
    truth_filename = replace_suffix(truth_path.stem, 'dqinit')
    truth_filename = truth_path.parent / (truth_filename + truth_path.suffix)

    diff = compare_asdf(dqinit_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_output_is_image_model(output_model):
    assert isinstance(output_model, rdm.datamodels.ImageModel)


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
    assert getattr(output_model.meta.cal_step, step_name) == "COMPLETE"


def test_has_a_wcs(output_model):
    assert output_model.meta.wcs is not None


@pytest.mark.soctests
def test_wcs_has_distortion_information(output_model):
    assert "v2v3" in output_model.meta.wcs.available_frames


def test_wcs_applies_distortion_correction(output_model):
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


@pytest.mark.parametrize(
    "arr_name", ("dq", "err", "var_poisson", "var_rnoise", "var_flat")
)
def test_array_exists(output_model, arr_name):
    assert hasattr(output_model, arr_name)


def test_has_exposure_time(output_model):
    assert "exposure_time" in output_model.meta.exposure


@pytest.mark.parametrize("meta_attribute", ("detector", "optical_element"))
def test_instrument_meta(output_model, meta_attribute):
    assert meta_attribute in output_model.meta.instrument


def test_wcs_has_bounding_box(output_model):
    assert len(output_model.meta.wcs.bounding_box) == 2
