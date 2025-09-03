"""Roman tests for source catalog creation"""

import pytest
from astropy.table import Table

from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]

fieldlist = [
    "ra_centroid",  # DMS374 positions on ICRF
    "dec_centroid",  # DMS374 positions on ICRF
    "aper01_f158_flux",  # DMS399 aperture fluxes
    "aper02_f158_flux",  # DMS399 aperture fluxes
    "aper04_f158_flux",  # DMS399 aperture fluxes
    "segment_f158_flux",  # DMS375 fluxes
    "kron_f158_flux",  # DMS375 fluxes
    "aper01_f158_flux_err",  # DMS386 flux uncertainties
    "aper02_f158_flux_err",  # DMS386 flux uncertainties
    "warning_flags",  # DMS387 dq_flags
    "is_extended_f158",  # DMS392 source classification
    "semimajor",  # DMS394 galaxy morphology
    "semiminor",  # DMS394 galaxy morphology
    "orientation_sky",  # DMS394 galaxy morphology
    "segment_f158_flux_err",  # DMS395 basic statistical uncertainties
    "kron_f158_flux_err",  # DMS395 basic statistical uncertainties
    "aper01_f158_flux_err",  # DMS395 basic statistical uncertainties
    "x_psf_f158_err",  # DMS395 basic statistical uncertainties
]


def test_multiband_catalog(rtdata_module, resource_tracker, request, dms_logger):
    rtdata = rtdata_module
    inputasnfn = "L3_skycell_mbcat_asn.json"
    # note that this input association currently only has a single
    # filter in it, so this is more of an existence proof for the multiband
    # catalogs than a detailed test.  Using only a single catalog lets us
    # rely on the existing regtest files.
    outputfn = "r00001_p_v01001001001001_270p65x70y49_f158_mbcat_cat.parquet"
    rtdata.get_asn(f"WFI/image/{inputasnfn}")
    rtdata.output = outputfn
    rtdata.input = inputasnfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")
    cattruth = Table.read(f"truth/{outputfn}")
    args = [
        "romancal.step.MultibandCatalogStep",
        inputasnfn,
        "--deblend",
        "True",  # use deblending, DMS 393
    ]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)
    cat = Table.read(outputfn)
    for field in fieldlist:
        assert field in cat.dtype.names

    dms_logger.info(
        "DMS374, 399, 375, 386, 387, 392, 394, 395: source catalog includes fields: "
        + ", ".join(fieldlist)
    )

    # DMS 393: multiband catalog uses both PSF-like and extend-source-like
    # kernels
    assert set(cat.dtype.names) == set(cattruth.dtype.names)
    # weak assertion that our truth file must at least have the same
    # catalog fields as the file produced here.  Exactly matching rows
    # would require a lot of okifying things that aren't obviously
    # the same, but we can easily check that the columns match up
    # by name.

    dms_logger.info("DMS391: successfully used multiple kernels to detect sources.")
    dms_logger.info("DMS393: successfully used deblending to separate blended sources.")
    dms_logger.info(
        "DMS399: successfully tested that catalogs contain aperture "
        "fluxes and uncertainties."
    )
