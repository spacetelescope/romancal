""" Roman tests for source catalog creation """

import asdf
import pytest

from romancal.source_catalog.source_catalog_step import SourceCatalogStep
from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]

input_filenames = [
    "r0099101001001001001_F158_visit_i2d.asdf",
    "r0000101001001001001_0001_wfi01_cal.asdf",
]

field_list = [
    "ra_centroid",  # DMS374 positions on ICRF
    "dec_centroid",  # DMS374 positions on ICRF
    "aper_total_flux",  # DMS375 fluxes
    "aper30_flux",  # DMS399 aperture fluxes
    "aper50_flux",  # DMS399 aperture fluxes
    "aper70_flux",  # DMS399 aperture fluxes
    "is_extended",  # DMS376 type of source
    "aper_total_flux_err",  # DMS386 flux uncertainties
    "flags",  # DMS387 dq_flags
]


@pytest.fixture(
    scope="module",
    params=input_filenames,
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
    field_list,
)
def test_has_field(fields, field):
    assert field in fields


def test_deblend_source_catalog(rtdata_module):
    rtdata = rtdata_module
    inputfn = input_filenames[1]  # use L2 image; could have picked either one
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn

    step = SourceCatalogStep()

    outputfn1 = inputfn.rsplit("_", 1)[0] + "_cat.asdf"
    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
        "--deblend",
        "False",
        "--output_file",
        outputfn1,
    ]
    RomanStep.from_cmdline(args)
    outputfn2 = inputfn.rsplit("_", 1)[0] + "_deblend_cat.asdf"
    rtdata.output = outputfn2
    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
        "--deblend",
        "True",
        "--output_file",
        outputfn2,
    ]
    RomanStep.from_cmdline(args)

    af1 = asdf.open(outputfn1)
    af2 = asdf.open(outputfn2)
    nsrc1 = len(af1["roman"]["source_catalog"])
    nsrc2 = len(af2["roman"]["source_catalog"])
    step.log.info(
        "DMS393: Deblended source catalog contains more sources "
        f"({nsrc2}) than undeblended ({nsrc1})?  "
        + ("PASS" if nsrc2 > nsrc1 else "FAIL")
    )
    assert nsrc2 > nsrc1


def test_kernel_detection(rtdata_module):
    rtdata = rtdata_module
    inputfn = input_filenames[1]  # use L2 image; could have picked either one
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    outputfn = inputfn.rsplit("_", 1)[0] + "_kernel10_cat.asdf"
    step = SourceCatalogStep()
    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
        "--kernel_fwhm",
        "10",
        "--output_file",
        outputfn,
    ]
    RomanStep.from_cmdline(args)
    af = asdf.open(outputfn)
    fields = af["roman"]["source_catalog"].dtype.names
    for field in field_list:
        assert field in field_list
    step.log.info("DMS391: Used alternative kernel to detect sources.  PASS")
