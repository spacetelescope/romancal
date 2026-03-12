"""Roman tests for source catalog creation"""

import numpy as np
import pytest
from astropy.table import Table
from roman_datamodels.datamodels import MultibandSegmentationMapModel

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
    "aper02_f129m_flux",  # DMS539 PSF-matched photometry
    "segment_f129m_flux",  # DMS539
    "kron_f129m_flux",  # DMS539
    "aper02_f213m_flux",  # DMS539 PSF-matched photometry
    "segment_f213m_flux",  # DMS539
    "kron_f213m_flux",  # DMS539
]


def test_multiband_catalog(rtdata_module, resource_tracker, request, dms_logger):
    rtdata = rtdata_module
    inputasnfn = "r00001_p_v01001001001001_270p65x70y49_asn.json"
    outputfn = "r00001_p_v01001001001001_270p65x70y49_cat.parquet"
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
        "--inject_sources",  # turn on source injection, DMS 396
        "True",
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

    # DMS 396: Ensure the segmentation image contains
    # both injected_sources and recovered_sources
    segm_mod = MultibandSegmentationMapModel(
        outputfn.replace("_cat.parquet", "_segm.asdf")
    )
    assert "injected_sources" in segm_mod
    assert "recovered_sources" in segm_mod

    dms_logger.info(
        "DMS396: segmentation image contains both "
        "injected_sources and recovered_sources."
    )

    # DMS 396: Ensure at least 50% of injected sourced are recovered.
    assert np.count_nonzero(segm_mod.recovered_sources["best_injected_index"] != -1) > (
        len(segm_mod.injected_sources) / 2
    )

    dms_logger.info("DMS396: successfully recovered over half of the injected sources.")

    # DMS 539/540: PSF matching quality check.
    # The aper02/aper08 flux ratio is a proxy for encircled energy and
    # should be similar across all matched bands (since they share an
    # effective PSF after matching).  Assert that matched ratios are much
    # closer to the F158 reference than the original unmatched ratios.
    bright = cat["segment_f158_flux"] > np.median(cat["segment_f158_flux"])
    ref_ratio = np.nanmedian(
        cat["aper02_f158_flux"][bright] / cat["aper08_f158_flux"][bright]
    )
    for band, matched_band in [("f129", "f129m"), ("f213", "f213m")]:
        unmatched_ratio = np.nanmedian(
            cat[f"aper02_{band}_flux"][bright] / cat[f"aper08_{band}_flux"][bright]
        )
        matched_ratio = np.nanmedian(
            cat[f"aper02_{matched_band}_flux"][bright]
            / cat[f"aper08_{matched_band}_flux"][bright]
        )
        dms_logger.info(
            f"aperture ratios: "
            f"{band}: ref={ref_ratio:.4f} matched={matched_ratio:.4f} "
            f"unmatched={unmatched_ratio:.4f}"
        )
        assert abs(matched_ratio - ref_ratio) < abs(unmatched_ratio - ref_ratio) / 5

    dms_logger.info(
        "DMS539: multiband catalog includes PSF-matched photometry for all bands."
    )
    dms_logger.info(
        "DMS540: PSF matching convolution kernels successfully constructed "
        "and used for each filter pair."
    )
