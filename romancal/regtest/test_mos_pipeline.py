"""Roman tests for the High Level Pipeline"""

import os

import pytest
import roman_datamodels as rdm

from romancal.pipeline.mosaic_pipeline import MosaicPipeline
from romancal.proj_match.proj_match import wcsinfo_to_wcs

from . import util
from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_mos(rtdata_module, resource_tracker):
    rtdata = rtdata_module

    rtdata.get_asn("WFI/image/L3_regtest_asn.json")

    # Test Pipeline
    output = "r0000101001001001001_f158_coadd.asdf"
    rtdata.output = output

    args = [
        "roman_mos",
        rtdata.input,
    ]
    with resource_tracker.track():
        MosaicPipeline.from_cmdline(args)

    rtdata.get_truth(f"truth/WFI/image/{output}")
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_mos):
    return run_mos.output


@pytest.fixture(scope="module")
def output_model(output_filename):
    with rdm.open(output_filename) as model:
        yield model


@pytest.fixture(scope="module")
def truth_filename(run_mos):
    return run_mos.truth


@pytest.fixture(scope="module")
def thumbnail_filename(output_filename):
    thumbnail_filename = output_filename.rsplit("_", 1)[0] + "_thumb.png"
    preview_cmd = f"stpreview to {output_filename} {thumbnail_filename} 256 256 roman"
    os.system(preview_cmd)  # noqa: S605
    return thumbnail_filename


@pytest.fixture(scope="module")
def preview_filename(output_filename):
    preview_filename = output_filename.rsplit("_", 1)[0] + "_preview.png"
    preview_cmd = f"stpreview to {output_filename} {preview_filename} 1080 1080 roman"
    os.system(preview_cmd)  # noqa: S605
    return preview_filename


def test_log_tracked_resources(log_tracked_resources, run_mos):
    log_tracked_resources()


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    # DMS356
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_thumbnail_exists(thumbnail_filename):
    # DMS356
    assert os.path.isfile(thumbnail_filename)


def test_preview_exists(preview_filename):
    # DMS356
    assert os.path.isfile(preview_filename)


@pytest.mark.parametrize("suffix", ("cat.parquet", "segm.asdf"))
def test_file_exists(output_filename, suffix):
    # DMS374 for catalog and segm
    expected_filename = output_filename.rsplit("_", 1)[0] + f"_{suffix}"
    assert os.path.isfile(expected_filename)


def test_output_is_mosaic(output_model):
    # DMS356
    assert isinstance(output_model, rdm.datamodels.MosaicModel)


@pytest.mark.parametrize(
    "step_name",
    (
        "skymatch",
        "outlier_detection",
        "resample",
    ),
)
def test_steps_ran(output_model, step_name):
    # DMS356
    # DMS400 for skymatch
    # DMS86 for outlier_detection and resample
    assert getattr(output_model.meta.cal_step, step_name) == "COMPLETE"


def test_added_background(output_model):
    # DMS400
    assert hasattr(output_model.meta.individual_image_meta, "background")


def test_added_background_level(output_model):
    # DMS400
    assert any(output_model.meta.individual_image_meta.background["level"] != 0)


def test_wcsinfo_wcs_roundtrip(output_model):
    """Test that the contents of wcsinfo reproduces the wcs"""
    wcs_from_wcsinfo = wcsinfo_to_wcs(output_model.meta.wcsinfo)

    ra_mad, dec_mad = util.comp_wcs_grids_arcs(output_model.meta.wcs, wcs_from_wcsinfo)
    assert (ra_mad + dec_mad) / 2.0 < 1.0e-5
