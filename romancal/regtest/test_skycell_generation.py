"""Roman tests for generating associatinos based on skycells"""

import json
import os

import pytest

from romancal.associations import skycell_asn

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]

EXPECTED_FILENAMES = [
    "r512_p_v01001001001_r274dp63x31y80_f158_asn.json",
    "r512_p_v01001001001_r274dp63x31y81_f158_asn.json",
    "r512_p_v01001001001_r274dp63x32y82_f158_asn.json",
    "r512_p_v01001001001_r274dp63x32y80_f158_asn.json",
    "r512_p_v01001001001_r274dp63x32y81_f158_asn.json",
    "r512_p_v01001001001_r274dp63x32y82_f158_asn.json",
]


@pytest.fixture(scope="module")
def run_skycell_asn(rtdata_module):
    rtdata = rtdata_module

    # This test should generate seven json files
    args = [
        "r0000101001001001001_0002_wfi01_cal.asdf",
        "r0000101001001001001_0002_wfi10_cal.asdf",
        "-o",
        "r512",
    ]
    rtdata.get_data("WFI/image/r0000101001001001001_0002_wfi01_cal.asdf")
    rtdata.get_data("WFI/image/r0000101001001001001_0002_wfi10_cal.asdf")

    skycell_asn._cli(args)
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
