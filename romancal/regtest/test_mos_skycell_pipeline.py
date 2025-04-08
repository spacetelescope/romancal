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

    # Test Pipeline
    rtdata.get_asn("WFI/image/L3_mosaic_asn.json")
    output = "r00001_p_v01001001001001_r274dp63x31y81_f158_coadd.asdf"
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


def test_log_tracked_resources(log_tracked_resources, run_mos):
    log_tracked_resources()


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_resample_ran(output_model):
    # DMS373 test output is resampled to a skycell
    # FIXME this doesn't test if the output is a skycell
    assert output_model.meta.cal_step.resample == "COMPLETE"


def test_location_name(output_model):
    # test that the location_name matches the skycell selected
    assert output_model.meta.basic.location_name == "r274dp63x31y81"


def test_wcsinfo_wcs_roundtrip(output_model):
    """Test that the contents of wcsinfo reproduces the wcs"""
    wcs_from_wcsinfo = wcsinfo_to_wcs(output_model.meta.wcsinfo)

    ra_mad, dec_mad = util.comp_wcs_grids_arcs(output_model.meta.wcs, wcs_from_wcsinfo)
    assert (ra_mad + dec_mad) / 2.0 < 1.0e-5
