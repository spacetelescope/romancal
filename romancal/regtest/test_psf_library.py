import asdf
import numpy as np
import pytest
from astropy import units as u
from roman_datamodels import datamodels as rdm

from romancal.step import SourceCatalogStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_psf_library(
    rtdata, ignore_asdf_paths, tmp_path, resource_tracker, dms_logger):

        input_data = "r0000101001001001001_0001_wfi01_f158_cal.asdf"

        rtdata.get_data(f"WFI/image/{input_data}")

        step = SourceCatalogStep()
        ref_file = step.get_reference_file(rtdata.input, "epsf")
        passmsg = "PASS" if ref_file not in"N/A"  else "FAIL"
        dms_logger.info(f"DMS531 {passmsg},  ePSF reference file {ref_file} exists.")
        passmsg = "PASS" if ref_file not in"N/A"  else "FAIL"
        dms_logger.info(f"DMS532 {passmsg},  ePSF reference file {ref_file} was retrieved.")


