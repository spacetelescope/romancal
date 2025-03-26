from pathlib import Path

import numpy as np
import pytest
import roman_datamodels as rdm
from numpy.testing import assert_allclose
from roman_datamodels.dqflags import pixel

from romancal.assign_wcs.assign_wcs_step import AssignWcsStep
from romancal.pipeline.exposure_pipeline import ExposurePipeline

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_elp(rtdata_module, resource_tracker):
    rtdata = rtdata_module

    input_data = "r0000201001001001001_0001_wfi01_uncal.asdf"
    rtdata.get_data(f"WFI/grism/{input_data}")
    rtdata.input = input_data

    # Test Pipeline
    output = "r0000201001001001001_0001_wfi01_cal.asdf"
    rtdata.output = output
    args = [
        "roman_elp",
        rtdata.input,
    ]
    with resource_tracker.track():
        ExposurePipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/grism/{output}")
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


def test_log_tracked_resources(log_tracked_resources, run_elp):
    log_tracked_resources()


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "step_name, status",
    (
        ("assign_wcs", "COMPLETE"),
        ("dark", "COMPLETE"),
        ("dq_init", "COMPLETE"),
        ("linearity", "COMPLETE"),
        ("ramp_fit", "COMPLETE"),
        ("saturation", "COMPLETE"),
        ("refpix", "COMPLETE"),
        # skipped
        ("flat_field", "SKIPPED"),
        ("photom", "SKIPPED"),
        ("source_catalog", "SKIPPED"),
        ("tweakreg", "SKIPPED"),
    ),
)
def test_step_status(output_model, step_name, status):
    # DMS90
    # DMS278
    # also DMS129 for assign_wcs
    assert getattr(output_model.meta.cal_step, step_name) == status


def test_jump_in_uneven_ramp(output_model):
    # DMS365 jump detection detected jumps in uneven ramp
    uneven = len({len(x) for x in output_model.meta.exposure.read_pattern}) > 1
    assert uneven & np.any(output_model.dq & pixel.JUMP_DET)


@pytest.mark.parametrize("arr_name", ("dq", "err", "var_poisson", "var_rnoise"))
def test_array_exists(output_model, arr_name):
    # DMS91
    assert hasattr(output_model, arr_name)


def test_has_exposure_time(output_model):
    # DMS88 total exposure time exists
    assert "exposure_time" in output_model.meta.exposure


def test_wcs_has_bounding_box(output_model):
    # DMS93 WCS tests
    assert len(output_model.meta.wcs.bounding_box) == 2


def test_repointed_matches_truth(
    repointed_filename_and_delta, rtdata, ignore_asdf_paths
):
    # DMS90
    repointed_filename, _ = repointed_filename_and_delta

    rtdata.get_truth(f"truth/WFI/grism/{Path(repointed_filename).name}")
    rtdata.output = repointed_filename
    diff = compare_asdf(repointed_filename, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_repointed_wcs_differs(repointed_filename_and_delta, output_model):
    # DMS93 WCS tests
    repointed_filename, delta = repointed_filename_and_delta
    orig_wcs = output_model.meta.wcs
    with rdm.open(repointed_filename) as repointed_model:
        assert_allclose(
            [orig_wcs(2048, 2048)[0] + delta[0], orig_wcs(2048, 2048)[1] + delta[1]],
            repointed_model.meta.wcs(2048, 2048),
            atol=1.0,
        )
