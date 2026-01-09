import pytest
import roman_datamodels as rdm

from romancal.pipeline.mosaic_pipeline import MosaicPipeline

from ..regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_mos(rtdata_module, resource_tracker):
    rtdata = rtdata_module

    rtdata.get_asn("WFI/image/compatibility/L3_regtest_asn.json")

    # Test Pipeline
    output = "r0000101001001001001_f158_coadd.asdf"
    rtdata.output = output

    args = [
        "roman_mos",
        rtdata.input,
    ]
    with resource_tracker.track():
        MosaicPipeline.from_cmdline(args)

    rtdata.get_truth(f"truth/WFI/image/compatibility/{output}")
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
