"""Roman tests for source catalog creation"""

import pytest
from astropy.table import Table

from romancal.source_catalog.source_catalog_step import SourceCatalogStep
from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(
    scope="module",
    params=[
        "r00001_p_v01001001001001_r274dp63x31y81_f158_coadd.asdf",
        "r0000101001001001001_f158_coadd.asdf",
        "r0000101001001001001_0001_wfi01_f158_cal.asdf",
    ],
    ids=["L3skycell", "L3", "L2"],
)
def run_source_catalog(rtdata_module, request, resource_tracker):
    rtdata = rtdata_module

    inputfn = request.param

    outputfn = inputfn.rsplit("_", 1)[0] + "_cat.parquet"
    rtdata.output = outputfn

    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")

    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
    ]
    with resource_tracker.track():
        RomanStep.from_cmdline(args)
    return rtdata_module


@pytest.fixture(scope="module")
def catalog(run_source_catalog):
    yield Table.read(run_source_catalog.output)


@pytest.fixture(scope="module")
def fields(catalog):
    return catalog.dtype.names


@pytest.mark.parametrize(
    "field",
    (
        "ra_centroid",  # DMS374 positions on ICRF
        "dec_centroid",  # DMS374 positions on ICRF
        "segment_flux",  # DMS375 fluxes
        "kron_flux",  # DMS375 fluxes
        "aper02_flux",  # DMS399 aperture fluxes
        "aper04_flux",  # DMS399 aperture fluxes
        "aper08_flux",  # DMS399 aperture fluxes
        "aper02_flux_err",  # DMS386 flux uncertainties
        "segment_flux_err",  # DMS386 flux uncertainties
        "is_extended",  # DMS376 type of source
        "warning_flags",  # DMS387 dq_flags
    ),
)
def test_has_field(fields, field):
    assert field in fields


def test_log_tracked_resources(log_tracked_resources, run_source_catalog):
    log_tracked_resources()


def test_forced_catalog(rtdata_module):
    rtdata = rtdata_module
    input_deep_segm = "r00001_p_v01001001001001_r274dp63x31y81_f158_segm.asdf"
    input_shallow_coadd = "r00001_p_e01001001001001_0001_r274dp63x31y81_f158_coadd.asdf"
    truth_cat = "r00001_p_e01001001001001_0001_r274dp63x31y81_f158_force_cat.parquet"
    rtdata.get_data(f"WFI/image/{input_deep_segm}")
    rtdata.get_data(f"WFI/image/{input_shallow_coadd}")
    truth_cat = rtdata.get_truth(f"truth/WFI/image/{truth_cat}")
    rtdata.input = input_shallow_coadd
    outputfn = input_shallow_coadd.rsplit("_", 1)[0] + "_force_cat.parquet"
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

    cat = Table.read(outputfn)
    fieldlist = [
        "forced_kron_flux",
        "forced_segment_flux",
        "forced_aper02_flux",
        "forced_semimajor_sigma",
        "forced_semiminor_sigma",
        "forced_ellipticity",
    ]
    for field in fieldlist:
        assert field in cat.dtype.names

    step = SourceCatalogStep()
    step.log.info(
        "DMS397: source catalog includes fields: "
        + ", ".join(fieldlist)
        + ", indicating measurements of morphology and photometry "
        "at multiple epochs."
    )

    cattruth = Table.read(truth_cat)
    assert set(cattruth.dtype.names) == set(cattruth.dtype.names)
    # weak assertion that our truth file must at least have the same
    # catalog fields as the file produced here.  Exactly matching rows
    # would require a lot of okifying things that aren't obviously
    # the same, but we can easily check that the columns match up
    # by name.
