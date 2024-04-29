""" Roman tests for source catalog creation """

import asdf
import pytest
from metrics_logger.decorators import metrics_logger

from romancal.source_catalog.source_catalog_step import SourceCatalogStep
from romancal.stpipe import RomanStep


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS374", "DMS375", "DMS376", "DMS386", "DMS387")
def test_catalog_l3(rtdata, ignore_asdf_paths):
    # DMS374: positions on ICRF
    # DMS375: fluxes
    # DMS376: type of source
    # DMS386: flux uncertainties
    # DMS387: DQ flags
    inputfn = "r0099101001001001001_F158_visit_0.900.0.50_178199.5_-0.5_i2d.asdf"
    outputfn = inputfn.replace("i2d", "cat")
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")

    outputfn = inputfn.replace("_i2d", "_cat")
    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
    ]
    RomanStep.from_cmdline(args)
    catalogfp = asdf.open(outputfn)
    catalog = catalogfp["roman"]["source_catalog"]
    step = SourceCatalogStep()

    fields = catalog.dtype.names
    has_pos = ("ra_centroid" in fields) and ("dec_centroid" in fields)
    step.log.info(
        "DMS374 MSG: L3 Catalog contains sky coordinates of sources? :"
        + passfail(has_pos)
    )
    assert has_pos
    has_flux = "aper_total_flux" in fields
    has_flux_err = "aper_total_flux_err" in fields
    has_type = "is_extended" in fields
    step.log.info("DMS375 MSG: L3 Catalog contains fluxes? :" + passfail(has_flux))
    assert has_flux
    step.log.info(
        "DMS376 MSG: L3 Catalog contains source classification? :" + passfail(has_type)
    )
    assert has_type
    step.log.info(
        "DMS386 MSG: L3 Catalog contains flux uncertainties? :" + passfail(has_flux_err)
    )
    assert has_flux_err

    # no compare_asdf on the catalogs
