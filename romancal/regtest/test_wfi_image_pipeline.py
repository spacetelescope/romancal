from pathlib import Path

import numpy as np
import pytest
import roman_datamodels as rdm
from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.dqflags import pixel

from romancal.assign_wcs.assign_wcs_step import AssignWcsStep
from romancal.pipeline.exposure_pipeline import ExposurePipeline

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = pytest.mark.bigdata


@pytest.fixture(scope="module")
def run_elp(rtdata_module):
    rtdata = rtdata_module

    input_data = "r0000101001001001001_0001_wfi01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000101001001001001_0001_wfi01_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
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


@pytest.fixture(scope="module")
def repointed_filename_and_delta(output_filename):
    delta = [10.0, 10.0]

    with rdm.open(output_filename) as model:
        # open from the file since we're going to modify the model
        del model.meta["wcs"]

        # Create new pointing for the model
        # RA & Dec are each shifted + 10 degrees, unless they are near
        # the upper limit, in which case they are shifted -10 degrees.
        if model.meta.wcsinfo.ra_ref >= 350.0:
            delta[0] *= -1.0
        if model.meta.wcsinfo.dec_ref >= 80.0:
            delta[1] *= -1.0

        model.meta.wcsinfo.ra_ref += delta[0]
        model.meta.wcsinfo.dec_ref += delta[1]

        # Create new wcs object for the new pointing
        model = AssignWcsStep.call(model)

        # save repointed model
        repointed_filename = output_filename.rsplit(".", 1)[0] + "_repoint.asdf"
        model.to_asdf(repointed_filename)

    return repointed_filename, delta


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
def test_jump_in_uneven_ramp(output_model):
    # DMS361 jump detection detected jumps in uneven ramp
    uneven = len({len(x) for x in output_model.meta.exposure.read_pattern}) > 1
    assert uneven & np.any(output_model.dq & pixel.JUMP_DET)


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


@pytest.mark.soctests
def test_repointed_matches_truth(
    repointed_filename_and_delta, rtdata, ignore_asdf_paths
):
    # DMS89 WCS tests
    repointed_filename, _ = repointed_filename_and_delta

    rtdata.get_truth(f"truth/WFI/image/{Path(repointed_filename).name}")
    rtdata.output = repointed_filename
    diff = compare_asdf(repointed_filename, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.soctests
def test_repointed_wcs_differs(repointed_filename_and_delta, output_model):
    # DMS89 WCS tests
    repointed_filename, delta = repointed_filename_and_delta
    orig_wcs = output_model.meta.wcs
    with rdm.open(repointed_filename) as repointed_model:
        assert_allclose(
            [orig_wcs(2048, 2048)[0] + delta[0], orig_wcs(2048, 2048)[1] + delta[1]],
            repointed_model.meta.wcs(2048, 2048),
            atol=1.0,
        )


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


@pytest.fixture(scope="module")
def run_all_saturated(rtdata_module):
    """
    Test ELP handling of an all saturated input model
    """
    rtdata = rtdata_module

    input_data = "r0000101001001001001_0001_wfi01_ALL_SATURATED_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000101001001001001_0001_wfi01_ALL_SATURATED_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
    ]
    ExposurePipeline.from_cmdline(args)

    # get truth file
    rtdata.get_truth(f"truth/WFI/image/{output}")
    return rtdata


@pytest.fixture(scope="module")
def all_saturated_model(run_all_saturated):
    with rdm.open(run_all_saturated.output) as model:
        yield model


def test_all_saturated_against_truth(run_all_saturated, ignore_asdf_paths):
    diff = compare_asdf(
        run_all_saturated.output, run_all_saturated.truth, **ignore_asdf_paths
    )
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "step_name, status",
    [
        ("dq_init", "COMPLETE"),
        ("saturation", "COMPLETE"),
        ("linearity", "SKIPPED"),
        ("dark", "SKIPPED"),
        ("ramp_fit", "SKIPPED"),
        ("assign_wcs", "SKIPPED"),
        ("flat_field", "SKIPPED"),
        ("photom", "SKIPPED"),
        ("source_catalog", "SKIPPED"),
        ("tweakreg", "SKIPPED"),
    ],
)
def test_all_saturated_status(all_saturated_model, step_name, status):
    """
    For an all saturated input the pipeline should skip all steps after saturation.
    """
    assert getattr(all_saturated_model.meta.cal_step, step_name) == status


def test_all_saturated_model_type(all_saturated_model):
    """
    For an all saturated input the output model should be an ImageModel.
    """
    assert isinstance(all_saturated_model, ImageModel)


@pytest.mark.parametrize("array_name", ["data", "err", "var_poisson", "var_rnoise"])
def test_all_saturated_zeroed(all_saturated_model, array_name):
    """
    For an all saturated input the output model should contain 0s for data and err arrays.
    """
    np.testing.assert_array_equal(getattr(all_saturated_model, array_name), 0)


def test_pipeline_suffix(rtdata, ignore_asdf_paths):
    """
    Tests passing suffix to the pipeline

    Note that this test mimics how the pipeline is run in OPS.

    Any changes to this test should be coordinated with OPS.
    """
    input_data = "r0000101001001001001_0001_wfi01_uncal.asdf"
    rtdata.get_data(f"WFI/image/{input_data}")

    output = "r0000101001001001001_0001_wfi01_star.asdf"
    rtdata.output = output

    args = [
        "roman_elp",
        rtdata.input,
        "--steps.tweakreg.skip=True",
        "--suffix=star",
    ]
    ExposurePipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

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
        assert model.meta.cal_step.source_catalog == "COMPLETE"
        assert model.meta.cal_step.tweakreg == "SKIPPED"
        assert model.meta.filename == output
