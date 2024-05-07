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


def check_catalog_fields(model, log, modeltype):
    fields = model.dtype.names
    has_pos = ("ra_centroid" in fields) and ("dec_centroid" in fields)
    log.info(
        f"DMS374 MSG: {modeltype} Catalog contains sky coordinates of sources? :"
        + passfail(has_pos)
    )

    has_flux = "aper_total_flux" in fields
    log.info(f"DMS375 MSG: {modeltype} Catalog contains fluxes? :" + passfail(has_flux))

    has_type = "is_extended" in fields
    log.info(
        f"DMS376 MSG: {modeltype} Catalog contains source classification? :"
        + passfail(has_type)
    )

    has_flux_err = "aper_total_flux_err" in fields
    log.info(
        f"DMS386 MSG: {modeltype} Catalog contains flux uncertainties? :"
        + passfail(has_flux_err)
    )

    has_flags = "flags" in fields
    log.info(
        f"DMS387 MSG: {modeltype} Catalog contains DQ flags? :" + passfail(has_flags)
    )

    return has_pos & has_flux & has_type & has_flux_err & has_flags


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS374", "DMS375", "DMS376", "DMS386", "DMS387")
def test_catalog_l3(rtdata, ignore_asdf_paths):
    # DMS374: positions on ICRF
    # DMS375: fluxes
    # DMS376: type of source
    # DMS386: flux uncertainties
    # DMS387: DQ flags
    inputfn = "r0099101001001001001_F158_visit_i2d.asdf"
    outputfn = inputfn.replace("_i2d", "_cat")
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")

    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
    ]
    RomanStep.from_cmdline(args)
    catalogfp = asdf.open(outputfn)
    catalog = catalogfp["roman"]["source_catalog"]
    step = SourceCatalogStep()
    assert check_catalog_fields(catalog, step.log, "L3")

    # no compare_asdf on the catalogs


@pytest.mark.bigdata
@pytest.mark.soctests
@metrics_logger("DMS374", "DMS375", "DMS376", "DMS386", "DMS387")
def test_catalog_l2(rtdata, ignore_asdf_paths):
    # DMS374: positions on ICRF
    # DMS375: fluxes
    # DMS376: type of source
    # DMS386: flux uncertainties
    # DMS387: DQ flags
    inputfn = "r0000101001001001001_01101_0001_WFI01_cal.asdf"
    outputfn = inputfn.replace("cal", "cat")
    rtdata.get_data(f"WFI/image/{inputfn}")
    rtdata.input = inputfn
    rtdata.get_truth(f"truth/WFI/image/{outputfn}")

    args = [
        "romancal.step.SourceCatalogStep",
        rtdata.input,
    ]
    RomanStep.from_cmdline(args)
    catalogfp = asdf.open(outputfn)
    catalog = catalogfp["roman"]["source_catalog"]
    step = SourceCatalogStep()

    assert check_catalog_fields(catalog, step.log, "L2")
    # no compare_asdf on the catalogs
