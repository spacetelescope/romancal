"""Roman tests for source catalog creation"""

import asdf
import pytest

from romancal.source_catalog.source_catalog_step import SourceCatalogStep
from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(
    scope="module",
    params=[
        "r0099101001001001001_r274dp63x31y81_prompt_F158_coadd.asdf",
        "r0099101001001001001_F158_visit_coadd.asdf",
        "r0000101001001001001_0001_wfi01_cal.asdf",
    ],
    ids=["L3", "L2", "L3skycell"],
)
def run_source_catalog(rtdata_module, request):
    rtdata = rtdata_module

    inputfn = request.param

    outputfn = inputfn.rsplit("_", 1)[0] + "_cat.asdf"
    rtdata.output = outputfn

    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")

    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
    ]
    RomanStep.from_cmdline(args)
    return rtdata_module


@pytest.fixture(scope="module")
def catalog(run_source_catalog):
    with asdf.open(run_source_catalog.output) as af:
        yield af["roman"]["source_catalog"]


@pytest.fixture(scope="module")
def fields(catalog):
    return catalog.dtype.names


@pytest.mark.parametrize(
    "field",
    (
        "ra_centroid",  # DMS374 positions on ICRF
        "dec_centroid",  # DMS374 positions on ICRF
        "aper30_flux",  # DMS399 aperture fluxes
        "aper50_flux",  # DMS399 aperture fluxes
        "aper70_flux",  # DMS399 aperture fluxes
        "aper_total_flux",  # DMS375 fluxes
        "is_extended",  # DMS376 type of source
        "aper_total_flux_err",  # DMS386 flux uncertainties
        "flags",  # DMS387 dq_flags
    ),
)
def test_has_field(fields, field):
    assert field in fields


def test_forced_catalog(rtdata_module):
    rtdata = rtdata_module
    input_deep_segm = "r0099101001001001001_r274dp63x31y81_prompt_F158_segm.asdf"
    input_shallow_coadd = (
        "r0099101001001001001_0001_r274dp63x31y81_prompt_F158_coadd.asdf"
    )
    truth_cat = "r0099101001001001001_0001_r274dp63x31y81_prompt_F158_force_cat.asdf"
    rtdata.get_data(f"WFI/image/{input_deep_segm}")
    rtdata.get_data(f"WFI/image/{input_shallow_coadd}")
    truth_cat = rtdata.get_truth(f"WFI/image/{truth_cat}")
    rtdata.input = input_shallow_coadd
    outputfn = input_shallow_coadd.rsplit("_", 1)[0] + "_force_cat.asdf"
    rtdata.output = outputfn

    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
        "--forced_segmentation",
        input_deep_segm,
        "--output_file",
        outputfn,
    ]
    RomanStep.from_cmdline(args)

    afcat = asdf.open(outputfn)
    fieldlist = [
        "forced_kron_flux",
        "forced_isophotal_flux",
        "forced_aper30_flux",
        "forced_semimajor_sigma",
        "forced_semiminor_sigma",
        "forced_ellipticity",
    ]
    for field in fieldlist:
        assert field in afcat["roman"]["source_catalog"].dtype.names

    step = SourceCatalogStep()
    step.log.info(
        "DMS397: source catalog includes fields: "
        + ", ".join(fieldlist)
        + ", indicating measurements of morphology and photometry "
        "at multiple epochs."
    )

    aftruth = asdf.open(truth_cat)
    assert set(aftruth["roman"]["source_catalog"].dtype.names) == set(
        afcat["roman"]["source_catalog"].dtype.names
    )
    # weak assertion that our truth file must at least have the same
    # catalog fields as the file produced here.  Exactly matching rows
    # would require a lot of okifying things that aren't obviously
    # the same, but we can easily check that the columns match up
    # by name.
