from pathlib import Path

import pytest

from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(scope="module")
def run_multiband_catalog(old_rtdata_module, request, resource_tracker):
    rtdata = old_rtdata_module

    # The filenames changed for B22
    _, build = rtdata._env.rsplit("_", maxsplit=1)
    if build <= "B21":
        input_fn = "L3_skycell_mbcat_asn.json"
        output_fn = "r00001_p_v01001001001001_270p65x70y49_f158_mbcat_cat.parquet"
    else:
        input_fn = "r00001_p_v01001001001001_270p65x70y49_asn.json"
        output_fn = "r00001_p_v01001001001001_270p65x70y49_cat.parquet"

    rtdata.get_asn(f"WFI/image/{input_fn}")

    rtdata.output = output_fn

    args = [
        "romancal.step.MultibandCatalogStep",
        rtdata.input,
        "--update_version=True",
    ]
    with resource_tracker.track():
        RomanStep.from_cmdline(args)
    return rtdata


def test_log_tracked_resources(log_tracked_resources, run_multiband_catalog):
    log_tracked_resources()


@pytest.fixture(scope="module")
def output_filename(run_multiband_catalog):
    return run_multiband_catalog.output


def test_output_exists(output_filename):
    assert Path(output_filename).exists
