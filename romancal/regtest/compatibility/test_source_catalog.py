from pathlib import Path

import pytest

from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(
    scope="module",
    params=[
        "r00001_p_v01001001001001_270p65x70y49_f158_coadd.asdf",
        "r0000101001001001001_f158_coadd.asdf",
        "r0000101001001001001_0001_wfi01_f158_cal.asdf",
    ],
    ids=["L3skycell", "L3", "L2"],
)
def run_source_catalog(rtdata_module, request, resource_tracker, old_build_path):
    rtdata = rtdata_module

    inputfn = request.param

    outputfn = inputfn.rsplit("_", 1)[0] + "_cat.parquet"
    rtdata.output = outputfn

    rtdata.get_data(f"{old_build_path}/WFI/image/compatibility/{inputfn}")
    rtdata.input = inputfn

    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
        "--update_version=True",
    ]
    with resource_tracker.track():
        RomanStep.from_cmdline(args)
    return rtdata_module


@pytest.fixture(scope="module")
def output_filename(run_source_catalog):
    return run_source_catalog.output


def test_log_tracked_resources(log_tracked_resources, run_source_catalog):
    log_tracked_resources()


def test_output_exists(output_filename):
    assert Path(output_filename).exists
