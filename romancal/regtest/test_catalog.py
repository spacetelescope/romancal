""" Roman tests for source catalog creation """

import asdf
import pytest

from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]


@pytest.fixture(
    scope="module",
    params=[
        "r0099101001001001001_F158_visit_i2d.asdf",
        "r0000101001001001001_0001_WFI01_cal.asdf",
    ],
    ids=["L3", "L2"],
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
        "aper_total_flux",  # DMS375 fluxes
        "is_extended",  # DMS376 type of source
        "aper_total_flux_err",  # DMS386 flux uncertainties
        "flags",  # DMS387 dq_flags
    ),
)
def test_has_field(fields, field):
    assert field in fields
