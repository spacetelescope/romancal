import numpy as np
from astropy.modeling.fitting import SplineSplrepFitter
from astropy.modeling.models import Spline1D
from roman_datamodels import datamodels
from roman_datamodels.datamodels import ImageModel, MosaicModel


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


def get_ee_spline(input_model, apcorr_file):
    """
    Create a spline fit to the encircled energy fraction vs radius data

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.ImageModel` or `~roman_datamodels.datamodels.MosaicModel`
        The input data model.
    """

    optical_element = input_model.meta.instrument.optical_element
    with datamodels.open(apcorr_file) as ee_ref:
        ee_fractions = getattr(ee_ref.data, optical_element).ee_fractions
        ee_radii = getattr(ee_ref.data, optical_element).ee_radii

        # Fit a spline model to the ee_fraction vs radius data so that we
        # can interpolate the values for arbitrary radii
        return SplineSplrepFitter()(Spline1D(), ee_radii, ee_fractions)


def copy_model_arrays(model):
    """
    Create a shallow copy of ImageModel or MosaicModel with data and err
    copied.

    This function creates a new model instance that shares the metadata
    with the input model but has independent copies of the data and err
    arrays. Other arrays (dq, weight) are shared references.

    Parameters
    ----------
    model : ImageModel or MosaicModel
        The input data model to copy.

    Returns
    -------
    copied_model : ImageModel or MosaicModel
        A new model with copied data and err arrays.

    Notes
    -----
    The metadata and dq/weight arrays are not copied because they are
    not modified in source catalog operations.
    """
    if isinstance(model, ImageModel):
        copied_model = ImageModel()
        copied_model.meta = model.meta
        copied_model.data = model.data.copy()
        copied_model.err = model.err.copy()
        copied_model.dq = model.dq
    elif isinstance(model, MosaicModel):
        copied_model = MosaicModel()
        copied_model.meta = model.meta
        copied_model.data = model.data.copy()
        copied_model.err = model.err.copy()
        copied_model.weight = model.weight
    else:
        raise TypeError("model must be an ImageModel or MosaicModel")

    return copied_model
