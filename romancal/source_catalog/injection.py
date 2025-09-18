"""
Module to inject sources into existing image / mosaic.
"""

import logging

from roman_datamodels.datamodels import ImageModel, MosaicModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

NOBJECTS = 300


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
    from romanisim.image import inject_sources_into_l2
    from romanisim.l3 import inject_sources_into_l3

    if isinstance(model, ImageModel):
        #  inject_sources_into_l2
        new_model = inject_sources_into_l2(model, si_cat)
    elif isinstance(model, MosaicModel):
        #  inject_sources_into_l3
        new_model = inject_sources_into_l3(model, si_cat)
    else:
        raise ValueError("The input model must be an ImageModel or MosaicModel.")

    return new_model


def make_cosmoslike_catalog(model, xpos, ypos):
    """
    Stuff

    Parameters
    ----------
    stuff

    Returns
    -------
    stuff
    """
    from romanisim.catalog import make_cosmos_galaxies

    # Pointing
    xcen, ycen = twcs.radecToxy(twcs.center.ra, twcs.center.dec, 'rad')

    # 75% of objects as galaxies
    # 25% of objects as point sources
    num_stars = int(len(xpos)/4)
    num_gals = len(xpos) - num_stars

    # Stars

    # Galaxies

