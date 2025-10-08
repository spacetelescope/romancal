"""
Module to inject sources into existing image / mosaic.
"""

import logging

import numpy as np
from astropy import table
from astropy import units as u
from roman_datamodels.datamodels import ImageModel, MosaicModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# One hour point source magnitude limit for each band
# https://roman.gsfc.nasa.gov/science/WFI_technical.html
HRPOINTMAGLIMIT = {
    "F062": 27.97,
    "F087": 27.63,
    "F106": 27.60,
    "F129": 27.60,
    "F158": 27.52,
    "F184": 26.95,
    "F213": 25.64,
    "F146": 28.01,
}

# One hour compact galaxy (0.3") source magnitude limit for each band
# https://roman.gsfc.nasa.gov/science/WFI_technical.html
HRGALMAGLIMIT = {
    "F062": 26.70,
    "F087": 26.38,
    "F106": 26.37,
    "F129": 26.37,
    "F158": 26.37,
    "F184": 25.95,
    "F213": 24.71,
    "F146": 26.84,
}


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


def make_cosmoslike_catalog(cen, xpos, ypos, exptime, filter="F146", seed=50, **kwargs):
    """
    Generate a catalog of cosmos galaxies and stars, with the following assumptions:
    - 75% of objects will be galaxies, 25% will be point sources
    - Galaxies will draw shapes and colors from the COSMOS catalog. Their position
    angles are drawn from a uniform distribution. Sersic indices to be drawn from
    between [0.8, 4.5]. Fluxes are scaled to the brightest color and vary in magnitude
    from 4 mag brighter to 1 mag fainter than the compact galaxy mag limit. Sizes scaled
    from magnitudes.
    - Stars will vary in magnitude from 4 mag brighter to 1 mag fainter than the point
    mag limit. They will have no color.

    Parameters
    ----------
    coord : astropy.coordinates.SkyCoord
        Location around which to generate sources.
    xpos, ypos : array_like (float)
        x, y positions of each source in objlist
    exptime : float
        Exposure time
    filter : str
        Name of filter (used to calculate size)
    seed : int
        Seed for random number generator

    Returns
    -------
    all_cat : astropy.Table
        Table for use with table_to_catalog to generate catalog for simulation.
    """
    from romanisim import bandpass, catalog

    # WFI bandpasses
    BANDPASSES = set(bandpass.galsim2roman_bandpass.values())

    # Set random source index for the catalog
    rng_numpy = np.random.default_rng(seed)
    ran_idx = rng_numpy.permutation(len(xpos))

    # 75% of objects as galaxies, 25% of objects as point sources
    num_stars = int(len(xpos) / 4)
    num_gals = len(xpos) - num_stars

    # Galaxies
    # RomanISim's make cosmos galaxies method will return cosmos like objects with
    # uniformly distributed position angles. Radius and area are set to return
    # a full sample (~346k galaxies)
    gal_cat = catalog.make_cosmos_galaxies(cen, radius=1.0, cat_area=(np.pi), **kwargs)

    # Trim to the required number of objects
    gal_cat = gal_cat[:num_gals]

    # Set magnitude spread
    mags = rng_numpy.uniform(low=-4.0, high=1.0, size=num_gals)

    # Brightest color flux
    gal_cat_data = (
        gal_cat.as_array(names=BANDPASSES)
        .view(dtype=float)
        .reshape((len(gal_cat), len(BANDPASSES)))
    )
    deepflux = gal_cat_data.max(axis=1)

    # Set bandpass fluxes
    for bp in BANDPASSES:
        # Normalize the mag limit to exptime
        gal_mag_limit = HRGALMAGLIMIT[bp] + (
            1.25 * np.log10((exptime * u.s).to(u.hour).value)
        )
        mag_tot = mags + gal_mag_limit

        # Set scaled fluxes
        gal_cat[bp] = (gal_cat[bp] * (10.0 ** (-(mag_tot) / 2.5)) / deepflux).astype(
            "f4"
        )

    # Sizes are drawn from a log-normal distribution of J-band magnitudes
    # The parameters below are derived from the distribution for sizes
    mu = (-0.1555 * ((-2.5 * np.log10(gal_cat["F129"])) - 17)) - 3.55
    sigma = 0.15
    radsize = rng_numpy.normal(mu, sigma)
    gal_cat["half_light_radius"] = (10**radsize) * 3600 * u.arcsec

    # Randomize concentrations
    gal_cat["n"] = rng_numpy.uniform(low=0.8, high=4.5, size=num_gals)

    # Stars
    # Create base table
    star_cat = table.Table()
    star_cat["ra"] = num_stars * [0]
    star_cat["dec"] = num_stars * [0]
    star_cat["type"] = num_stars * ["PSF"]
    star_cat["n"] = num_stars * [-1]
    star_cat["half_light_radius"] = num_stars * [0]
    star_cat["pa"] = num_stars * [0]
    star_cat["ba"] = num_stars * [1]

    # Set magnitude spread
    mags = rng_numpy.uniform(low=-4.0, high=1.0, size=num_stars)

    # Set bandpass fluxes
    for bp in BANDPASSES:
        point_mag_limit = HRPOINTMAGLIMIT[bp] + (
            1.25 * np.log10((exptime * u.s).to(u.hour).value)
        )
        star_cat[bp] = (10.0 ** (-(mags + point_mag_limit) / 2.5)).astype("f4")

    # Combine the objects
    all_cat = table.vstack([gal_cat, star_cat])

    # Set the positions to randomly selected objects
    all_cat["ra"] = np.array(xpos)[ran_idx].tolist()
    all_cat["dec"] = np.array(ypos)[ran_idx].tolist()

    return all_cat
