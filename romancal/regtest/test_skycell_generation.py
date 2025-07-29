"""Roman tests for generating associatinos based on skycells"""

import json
import os

import pytest
import roman_datamodels.datamodels as rdm

from romancal.associations import skycell_asn
from romancal.pipeline.mosaic_pipeline import MosaicPipeline

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]

EXPECTED_FILENAMES = [
    "r00001_p_v01001001001001_270p65x69y48_f158_asn.json",
    "r00001_p_v01001001001001_270p65x70y48_f158_asn.json",
    "r00001_p_v01001001001001_270p65x71y48_f158_asn.json",
    "r00001_p_v01001001001001_270p65x69y49_f158_asn.json",
    "r00001_p_v01001001001001_270p65x70y49_f158_asn.json",
    "r00001_p_v01001001001001_270p65x71y49_f158_asn.json",
    "r00001_p_v01001001001001_270p65x69y50_f158_asn.json",
    "r00001_p_v01001001001001_270p65x70y50_f158_asn.json",
    "r00001_p_v01001001001001_270p65x71y50_f158_asn.json",
    "r00001_p_v01001001001001_270p65x69y51_f158_asn.json",
    "r00001_p_v01001001001001_270p65x70y51_f158_asn.json",
    "r00001_p_v01001001001001_270p65x71y51_f158_asn.json",
    "r00001_p_v01001001001001_270p65x69y52_f158_asn.json",
    "r00001_p_v01001001001001_270p65x70y52_f158_asn.json",
    "r00001_p_v01001001001001_270p65x71y52_f158_asn.json",
]


@pytest.fixture(scope="module")
def run_skycell_asn(rtdata_module):
    rtdata = rtdata_module

    # This test should generate seven json files
    args = [
        "r0000101001001001001_0002_wfi01_f158_cal.asdf",
        "r0000101001001001001_0002_wfi10_f158_cal.asdf",
        "-o",
        "r00001",
    ]
    rtdata.get_data("WFI/image/r0000101001001001001_0002_wfi01_f158_cal.asdf")
    rtdata.get_data("WFI/image/r0000101001001001001_0002_wfi10_f158_cal.asdf")

    skycell_asn._cli(args)
    return rtdata


@pytest.fixture(scope="module")
def mosaic_pipeline_on_skycell_asn(run_skycell_asn):
    rtdata = run_skycell_asn
    rtdata.input = EXPECTED_FILENAMES[0]
    rtdata.output = f"{rtdata.input.rsplit('_', maxsplit=1)[0]}_coadd.asdf"

    # just run the first association (not all)
    args = [
        "roman_mos",
        rtdata.input,
    ]

    # we don't setup or fetch a truth file here as the aim of this
    # test is to check the output is resampled onto a skycell
    MosaicPipeline.from_cmdline(args)
    return rtdata


@pytest.mark.parametrize("expected_filename", EXPECTED_FILENAMES)
def test_file_exists(run_skycell_asn, expected_filename):
    """Test that the expected json files were generated"""
    assert os.path.isfile(expected_filename)


@pytest.mark.parametrize("expected_filename", EXPECTED_FILENAMES)
def test_files_contain_wcsinfo(run_skycell_asn, expected_filename):
    with open(expected_filename) as f:
        asn = json.load(f)
    assert "skycell_wcs_info" in asn


def test_mosaic_output_is_skycell(mosaic_pipeline_on_skycell_asn):
    """Check that the mos output for the generated association is on a skycell"""
    rtdata = mosaic_pipeline_on_skycell_asn
    filename = rtdata.output
    with rdm.open(filename) as model:
        assert model.meta.wcsinfo.skycell_name == "270p65x48y69"
