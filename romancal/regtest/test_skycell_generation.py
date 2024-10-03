""" Roman tests for generating associatinos based on skycells"""

import os

import pytest

from romancal.associations import skycell_asn

def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


@pytest.mark.bigdata
def test_skycell_asn_generation(rtdata, ignore_asdf_paths):
    """Test for the generation of associations based on skycells"""

    # This test should generate seven json files
    args = ['r0000101001001001001_01101_0002_WFI01_cal.asdf', \
    'r0000101001001001001_01101_0002_WFI10_cal.asdf', '-o', 'r512']
    #rtdata.get_asn("WFI/image/L3_regtest_asn.json")
    rtdata.get_data("WFI/image/r0000101001001001001_01101_0002_WFI01_cal.asdf")
    rtdata.get_data("WFI/image/r0000101001001001001_01101_0002_WFI10_cal.asdf")

    skycell_asn.Main(args)

    output_files = ['r512_r274dp63x31y80_visit_F158_prompt_i2d_asn.json',
                    'r512_r274dp63x31y81_visit_F158_prompt_i2d_asn.json',
                    'r512_r274dp63x31y82_visit_F158_prompt_i2d_asn.json',
                    'r512_r274dp63x32y80_visit_F158_prompt_i2d_asn.json',
                    'r512_r274dp63x32y81_visit_F158_prompt_i2d_asn.json',
                    'r512_r274dp63x33y80_visit_F158_prompt_i2d_asn.json',
                    'r512_r274dp63x33y81_visit_F158_prompt_i2d_asn.json']
    # Test that the json files exist
    for file in output_files:
        skycell_asn.logger.info(
            "Check that the json file exists   "
            + passfail(os.path.isfile(file))
        )
