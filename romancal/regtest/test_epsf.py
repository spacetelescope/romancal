"""Regression test for the ePSF reference file script."""

import numpy as np
import pytest
import roman_datamodels as rdm

from romancal.source_catalog import _epsf
from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = pytest.mark.bigdata

INPUT = "r0000101001001001001_0001_wfi01_f158_cal.asdf"
CATALOG = "r0000101001001001001_0001_wfi01_f158_cat.parquet"
OUTPUT = "r0000101001001001001_0001_wfi01_f158_epsf.asdf"


@pytest.fixture(scope="module")
def run_epsf(rtdata_module, resource_tracker):
    """Make a catalog for the L2 image, then run roman_epsf_from_model."""
    rtdata = rtdata_module
    rtdata.get_data(f"WFI/image/{INPUT}")
    rtdata.input = INPUT
    rtdata.output = OUTPUT
    rtdata.get_truth(f"truth/WFI/image/{OUTPUT}")

    # Make the catalog here so that it has the columns pick_stars uses.
    RomanStep.from_cmdline(["romancal.step.SourceCatalogStep", INPUT])
    with resource_tracker.track():
        # The output name defaults to OUTPUT.
        _epsf.cli([INPUT, CATALOG])
    with rdm.open(rtdata.output) as model, rdm.open(rtdata.truth) as truth:
        yield model, truth


def test_epsf_validates(run_epsf, dms_logger):
    """The output is a valid ePSF reference file for the input image."""
    epsf_model, _ = run_epsf
    # DMS 530 is to construct empirical PSFs. The script makes an ePSF
    # reference file for the detector and filter of the input image.
    epsf_model.validate()
    matches = (
        epsf_model.meta.instrument.detector == "WFI01"
        and epsf_model.meta.instrument.optical_element == "F158"
    )
    passmsg = "PASS" if matches else "FAIL"
    dms_logger.info(
        f"DMS530 {passmsg}, empirical ePSF reference file is valid and "
        "matches the input detector and filter."
    )
    assert matches


def test_epsf_values(run_epsf, dms_logger):
    """The ePSF is finite, centered, and sums to oversample**2."""
    epsf_model, _ = run_epsf
    # DMS 530 is to construct empirical PSFs. Check that the empirical
    # ePSF is a sensible PSF.
    oversample = epsf_model.meta.oversample
    psf = epsf_model.psf[0, 0, 0]
    peak = np.unravel_index(np.argmax(psf), psf.shape)
    center = (np.array(psf.shape) - 1) / 2
    checks = {
        "finite": np.all(np.isfinite(psf)),
        "normalized": np.isclose(psf.sum(), oversample**2, rtol=1e-2),
        "centered": np.all(np.array(peak) == center),
    }
    passed = all(checks.values())
    passmsg = "PASS" if passed else "FAIL"
    dms_logger.info(
        f"DMS530 {passmsg}, empirical ePSF is finite, normalized, and centered."
    )
    assert passed, checks


def test_epsf_matches_truth(run_epsf):
    """The ePSF has the same shape and oversampling as the truth file.

    The values are not compared; they change with small changes in the
    star selection.
    """
    epsf_model, truth = run_epsf
    assert epsf_model.psf.shape == truth.psf.shape
    assert epsf_model.meta.oversample == truth.meta.oversample
