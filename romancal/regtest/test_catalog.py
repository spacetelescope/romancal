""" Roman tests for source catalog creation """

import asdf
import pytest

from romancal.stpipe import RomanStep

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]

input_filenames = [
    "r0099101001001001001_F158_visit_i2d.asdf",
    "r0000101001001001001_0001_wfi01_cal.asdf",
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
    (
        "ra_centroid",  # DMS374 positions on ICRF
        "dec_centroid",  # DMS374 positions on ICRF
        "aper_total_flux",  # DMS375 fluxes
        "aper30_flux",  # DMS399 aperture fluxes
        "aper50_flux",  # DMS399 aperture fluxes
        "aper70_flux",  # DMS399 aperture fluxes
        "is_extended",  # DMS376 type of source
        "aper_total_flux_err",  # DMS386 flux uncertainties
        "flags",  # DMS387 dq_flags
    ),
)
def test_has_field(fields, field):
    assert field in fields


def test_deblend_source_catalog(rtdata_module):
    rtdata = rtdata_module
    inputfn = input_filenames[1]
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn

    outputfn1 = inputfn.rsplit("_", 1)[0] + "_cat.asdf"
    args = ["romancal.step.SourceCatalogStep", rtdata.input,
            "--deblend", "False", "--output_file", outputfn1]
    RomanStep.from_cmdline(args)
    outputfn2 = inputfn.rsplit("_", 1)[0] + "_deblend_cat.asdf"
    rtdata.output = outputfn2
    args = ["romancal.step.SourceCatalogStep", rtdata.input,
            "--deblend", "True", "--output_file", outputfn2]
    RomanStep.from_cmdline(args)

    af1 = asdf.open(outputfn1)
    af2 = asdf.open(outputfn2)
    nsrc1 = len(af1['roman']['source_catalog'])
    nsrc2 = len(af2['roman']['source_catalog'])
    from romancal.source_catalog.source_catalog_step import SourceCatalogStep
    step = SourceCatalogStep()
    step.log.info(
        "DMS393: Deblended source catalog contains more sources "
        f"({nsrc2}) than undeblended ({nsrc1})?  "
        + ('PASS' if nsrc2 > nsrc1 else 'FAIL'))
    print('hello!')
    # can't figure out where these log messages are going??
    assert nsrc2 > nsrc1
    # what's the nicer way to do this with fixtures?
    

def test_kernel_detection(rtdata_module):
    rtdata = rtdata_module
    inputfn = input_filenames[1]
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    outputfn = inputfn.rsplit("_", 1)[0] + "_kernel10_cat.asdf"
    args = ["romancal.step.SourceCatalogStep", rtdata.input,
            "--kernel_fwhm", "20", "--output_file", outputfn]
    RomanStep.from_cmdline(args)
    af = asdf.open(outputfn)
    # test_has_fields on this catalog with the large kernel?
