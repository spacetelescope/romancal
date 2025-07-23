import numpy as np


def copy_mosaic_meta(model, cat_model):
    # TODO some junk values here
    cat_model.meta.prd_version = "8.8.8"
    cat_model.meta.sdf_software_version = "7.7.7"
    cat_model.meta.basic.time_first_mjd = model.meta.coadd_info.time_first.mjd
    cat_model.meta.basic.time_last_mjd = model.meta.coadd_info.time_last.mjd
    cat_model.meta.basic.time_mean_mjd = model.meta.coadd_info.time_mean.mjd
    cat_model.meta.basic.max_exposure_time = model.meta.coadd_info.get(
        "max_exposure_time", np.nan
    )
    cat_model.meta.basic.mean_exposure_time = model.meta.coadd_info.exposure_time
    cat_model.meta.basic.visit = model.meta.observation.visit
    cat_model.meta.basic.segment = model.meta.observation.segment
    cat_model.meta.basic["pass"] = model.meta.observation["pass"]
    cat_model.meta.basic.program = model.meta.observation.program
    cat_model.meta.basic.survey = "?"

    # TODO handle optical element with multiple values
    cat_model.meta.basic.optical_element = model.meta.instrument.optical_element

    cat_model.meta.basic.instrument = "WFI"
    cat_model.meta.basic.location_name = model.meta.wcsinfo.skycell_name
    cat_model.meta.basic.product_type = model.meta.product_type

    # TODO can we fill these in?
    cat_model.meta.photometry.conversion_megajanskys = None
    cat_model.meta.photometry.conversion_megajanskys_uncertainty = None
    cat_model.meta.photometry.pixel_area = None
    return cat_model
