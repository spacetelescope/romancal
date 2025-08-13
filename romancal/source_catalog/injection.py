"""
Module to inject sources into existing image / mosaic.
"""

import logging
from romanisim.image import inject_sources_into_l2
from romanisim.l3 import inject_sources_into_l3
from roman_datamodels.datamodels import ImageModel, MosaicModel


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def inject_sources(model, si_cat):
    """
    Convolve the background-subtracted model image with a Gaussian
    kernel.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        Model into which to inject sources.
    si_cat: astropy.table.Table
        Catalog of sources to inject into image.

    Returns
    -------
    new_model : `ImageModel` or `MosaicModel`
        Input model with added sources.
    """
    if isinstance(model, ImageModel):
        #  inject_sources_into_l2
        new_model = inject_sources_into_l2(model, si_cat)
    elif isinstance(model, MosaicModel):
        #  inject_sources_into_l3
        new_model = model.copy()
        _ = inject_sources_into_l3(new_model, si_cat)
    else:
        raise ValueError("The input model must be an ImageModel or MosaicModel.")

    return new_model
