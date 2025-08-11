"""Tests for the Roman PSF Library"""

import numpy as np
import pytest
from roman_datamodels import datamodels as rdm

from romancal.source_catalog import psf
from romancal.step import SourceCatalogStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(
    scope="module",
    params=[
        "r0000101001001001001_0001_wfi01_f158_cal.asdf",
    ],
    ids=["input_file"],
)
def get_input_file(rtdata_module, request, resource_tracker):
    """Get the input file for testing"""

    rtdata = rtdata_module

    input_file = request.param

    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file
    return rtdata


@pytest.mark.bigdata
def test_psf_library_reffile(get_input_file, dms_logger):
    """Test the retrieval of the PSF reference file"""

    rtdata = get_input_file

    # DMS 531 is to check that a PSF library has been provided and
    # DMS 532 is that we have access to the PSF library. By retrieving a
    # PSF reference file we show that the library exists and we have access.
    step = SourceCatalogStep()
    ref_file = step.get_reference_file(rtdata.input, "epsf")
    has_ref_file = ref_file != "N/A"
    assert has_ref_file
    passmsg = "PASS" if has_ref_file else "FAIL"
    dms_logger.info(f"DMS531 {passmsg},  ePSF reference file {ref_file} exists.")
    dms_logger.info(f"DMS532 {passmsg},  ePSF reference file {ref_file} was retrieved.")


@pytest.mark.bigdata
def test_psf_library_crdsfile(get_input_file, dms_logger):
    """Test that the PSF reference file matches the observation"""
    # DMS 535 identify an appropriate PSF for WFI source
    # check that the detector and optical element for the psf ref file
    # matches the data file

    rtdata = get_input_file
    input_data = rtdata.input

    step = SourceCatalogStep()
    ref_file = step.get_reference_file(input_data, "epsf")
    with rdm.open(ref_file) as ref_data:
        with rdm.open(rtdata.input) as wfi_data:
            psf_selection_match = (
                wfi_data.meta.instrument.optical_element
                == ref_data.meta.instrument.optical_element
            ) and (
                wfi_data.meta.instrument.detector == ref_data.meta.instrument.detector
            )

            assert psf_selection_match
            passmsg = "PASS" if psf_selection_match else "FAIL"
            dms_logger.info(
                f"DMS535 {passmsg},  ePSF reference file selection matches input data file"
            )


@pytest.mark.bigdata
def test_psf_library_psfinterp(get_input_file, dms_logger):
    """Test that the interpolation is occurring"""
    # DMS 536 interpolating empirical ePSFs given the observed-source position
    # Create two psf's one at the center and one off-center to show we can interpolate the PSF
    # and that the two are not the same

    rtdata = get_input_file
    input_data = rtdata.input

    step = SourceCatalogStep()
    ref_file = step.get_reference_file(input_data, "epsf")
    with rdm.open(ref_file) as ref_data:
        psf_model = psf.get_gridded_psf_model(ref_data)
        center = 2044
        szo2 = 10  # psf stamp axis lenght/2
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
        assert psf_diff
        passmsg = "PASS" if psf_diff else "FAIL"
        dms_logger.info(
            f"DMS536 {passmsg},  The interpolated ePSF for two positions are not the same"
        )
