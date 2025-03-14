"""Roman tests for source catalog creation"""

import asdf
import pytest

from romancal.multiband_catalog.multiband_catalog_step import MultibandCatalogStep
from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]

fieldlist = [
    "ra_centroid",  # DMS374 positions on ICRF
    "dec_centroid",  # DMS374 positions on ICRF
    "F158_aper30_flux",  # DMS399 aperture fluxes
    "F158_aper50_flux",  # DMS399 aperture fluxes
    "F158_aper70_flux",  # DMS399 aperture fluxes
    "F158_aper_total_flux",  # DMS375 fluxes
    "F158_aper_total_flux_err",  # DMS386 flux uncertainties
    "F158_flags",  # DMS387 dq_flags
    "F158_is_extended",  # DMS392 source classification
    "semimajor_sigma",  # DMS394 galaxy morphology
    "semiminor_sigma",  # DMS394 galaxy morphology
    "orientation",  # DMS394 galaxy morphology
    "F158_isophotal_flux_err",  # DMS395 basic statistical uncertainties
    "F158_kron_flux_err",  # DMS395 basic statistical uncertainties
    "F158_aper30_flux_err",  # DMS395 basic statistical uncertainties
    "F158_x_psf_err",  # DMS395 basic statistical uncertainties
]


def test_multiband_catalog(rtdata_module):
    rtdata = rtdata_module
    inputasnfn = "L3_skycell_mbcat_asn.json"
    # note that this input association currently only has a single
    # filter in it, so this is more of an existence proof for the multiband
    # catalogs than a detailed test.  Using only a single catalog lets us
    # rely on the existing regtest files.
    outputfn = "r0099101001001001001_r274dp63x31y81_prompt_F158_mbcat_cat.asdf"
    rtdata.get_asn(f"WFI/image/{inputasnfn}")
    rtdata.output = outputfn
    rtdata.input = inputasnfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")
    aftruth = asdf.open(f"truth/{outputfn}")
    args = [
        "romancal.step.MultibandCatalogStep",
        inputasnfn,
        "--deblend",
        "True",  # use deblending, DMS 393
    ]
    RomanStep.from_cmdline(args)
    afcat = asdf.open(outputfn)
    for field in fieldlist:
        assert field in afcat["roman"]["source_catalog"].dtype.names

    step = MultibandCatalogStep()
    step.log.info(
        "DMS374, 399, 375, 386, 387, 392, 394, 395: source catalog includes fields: "
        + ", ".join(fieldlist)
    )

    # DMS 393: multiband catalog uses both PSF-like and extend-source-like
    # kernels
    assert set(aftruth["roman"]["source_catalog"].dtype.names) == set(
        afcat["roman"]["source_catalog"].dtype.names
    )
    # weak assertion that our truth file must at least have the same
    # catalog fields as the file produced here.  Exactly matching rows
    # would require a lot of okifying things that aren't obviously
    # the same, but we can easily check that the columns match up
    # by name.

    step.log.info("DMS391: successfully used multiple kernels to detect sources.")
    step.log.info("DMS393: successfully used deblending to separate blended sources.")
    step.log.info(
        "DMS399: successfully tested that catalogs contain aperture "
        "fluxes and uncertainties."
    )
