import numpy as np
import pytest

from roman_datamodels import datamodels as rdm

import romancal.source_catalog.psf as psf
from romancal.step import SourceCatalogStep


@pytest.mark.bigdata
def test_psf_library(
    rtdata, ignore_asdf_paths, tmp_path, resource_tracker, dms_logger):

        input_data = "r0000101001001001001_0001_wfi01_f158_cal.asdf"

        rtdata.get_data(f"WFI/image/{input_data}")
        rtdata.input = input_data

        step = SourceCatalogStep()
        # DMS 531 is to check that a PSF library has been provided and
        # DMS 532 is that we have access to the PSF library. By retrieving a
        # PSF reference file we show that the library exists and we have access. 
        ref_file = step.get_reference_file(rtdata.input, "epsf")
        has_ref_file = ref_file != "N/A"
        assert has_ref_file
        passmsg = "PASS" if has_ref_file  else "FAIL"
        dms_logger.info(f"DMS531 {passmsg},  ePSF reference file {ref_file} exists.")
        passmsg = "PASS" if has_ref_file  else "FAIL"
        dms_logger.info(f"DMS532 {passmsg},  ePSF reference file {ref_file} was retrieved.")

        # DMS 535 identify an appropriate PSF for WFI source
        # check that the detector and optical element for the psf ref file
        # matches the data file
        with rdm.open(ref_file) as ref_data:
                with rdm.open(rtdata.input) as wfi_data:
                        psf_selection_match = ((wfi_data.meta.instrument.optical_element == ref_data.meta.instrument.optical_element) and
                                              (wfi_data.meta.instrument.detector == ref_data.meta.instrument.detector))
                        assert psf_selection_match
                        passmsg = "PASS" if psf_selection_match else "FAIL"
                        dms_logger.info(f"DMS535 {passmsg},  ePSF reference file selection matches input data file")                             

                        # DMS 536 interpolating empirical ePSFs given the observed-source position
                        # Create two psf's one at the center and one off-center to show we can interpolate the PSF
                        # and that the two are not the same
                        psf_model = psf.get_gridded_psf_model(ref_data)
                        center = 2044
                        szo2 = 10 # psf stamp axis lenght/2
                        oversample = 11
                        npix = szo2 * 2 * oversample + 1
                        pts = np.linspace(-szo2 + center, szo2 + center, npix)
                        xx, yy = np.meshgrid(pts, pts)
                        psf_center = psf_model.evaluate(xx, yy, 1, center, center)
                        off_center = center + 99
                        pts = np.linspace(-szo2 + off_center, szo2 + off_center, npix)
                        xx, yy = np.meshgrid(pts, pts)
                        psf_offcenter = psf_model.evaluate(xx, yy, 1, off_center, off_center)
                        # show that these two are different
                        diff = psf_offcenter - psf_center
                        # check to make sure the difference is not zero over the array
                        psf_diff = diff.any()
                        passmsg = "PASS" if psf_diff else "FAIL"
                        dms_logger.info(f"DMS536 {passmsg},  The interpolated ePSF for two positions are not the same") 
