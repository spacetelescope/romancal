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
def render_psfs(rtdata_module, request, resource_tracker):
    """Get the input file for testing"""

    rtdata = rtdata_module

    input_file = request.param

    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file
    rtdata.output = "psf_render.asdf"
    rtdata.get_truth(f"truth/WFI/image/psf_render.asdf")

    step = SourceCatalogStep()
    ref_file = step.get_reference_file(rtdata.input, "epsf")
    grid_nominal = psf.get_gridded_psf_model(ref_data, 0, 1)
    grid_defocus = psf.get_gridded_psf_model(ref_data, 1, 1)
    grid_red = psf.get_gridded_psf_model(ref_data, 0, 2)
    out = dict()
    pixcen = 2044
    out['stamp_center'] = psf.render_stamp(pixcen, pixcen, grid, 19)
    out['stamp_corner'] = psf.render_stamp(0, 0, grid, 19)
    out['stamp_red'] = psf.render_stamp(pixcen, pixcen, grid_red, 19)
    out['stamp_defocus'] = psf.render_stamp(pixcen, pixcen, grid_defocus, 19)
    asdf.dump(out, open(rtdata.output, 'wb'))
    return rtdata, out


@pytest.mark.bigdata
def test_psf_library_reffile(render_psfs, dms_logger):
    """Test the retrieval of the PSF reference file"""

    rtdata, stamps = render_psfs

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
def test_psf_library_crdsfile(render_psfs, dms_logger):
    """Test that the PSF reference file matches the observation"""
    # DMS 535 identify an appropriate PSF for WFI source
    # check that the detector and optical element for the psf ref file
    # matches the data file

    rtdata, stamps = render_psfs
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
def test_psf_library_psfinterp(render_psfs, dms_logger):
    """Test that the interpolation is occurring"""
    # DMS 536 interpolating empirical ePSFs given the observed-source position
    # Create two psf's one at the center and one off-center to show we can interpolate the PSF
    # and that the two are not the same

    _, stamps = render_psfs

    diff = stamps['stamp_center'] - stamps['stamp_corner']
    # check to make sure the difference is not zero over the array
    psf_diff = diff.any()
    assert psf_diff
    passmsg = "PASS" if psf_diff else "FAIL"
    dms_logger.info(
        f"DMS536 {passmsg},  The interpolated ePSF for two positions are not the same"
    )

@pytest.mark.bigdata
def test_psf_library_variability(render_psfs, dms_logger):
    """Test that PSF variation with color is tracked"""
    # DMS 533 handling variable PSFs
    # Create two PSFs, one for a red source, one for a blue source, show that they're
    # different

    _, stamps = render_psfs

    diff = stamps['stamp_red'] - stamps['stamp_center']
    # check to make sure the difference is not zero over the array
    psf_diff = diff.any()
    assert psf_diff
