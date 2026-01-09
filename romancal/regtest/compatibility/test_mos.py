from pathlib import Path

import pytest

from romancal.pipeline.mosaic_pipeline import MosaicPipeline

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_mos(old_rtdata_module, resource_tracker):
    rtdata = old_rtdata_module

    rtdata.get_asn("WFI/image/L3_regtest_asn.json")

    # Test Pipeline
    output = "r0000101001001001001_f158_coadd.asdf"
    rtdata.output = output

    args = [
        "roman_mos",
        rtdata.input,
        "--update_version=True",
    ]
    with resource_tracker.track():
        MosaicPipeline.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_mos):
    return run_mos.output


def test_log_tracked_resources(log_tracked_resources, run_mos):
    log_tracked_resources()


def test_output_exists(output_filename):
    assert Path(output_filename).exists()
