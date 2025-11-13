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


def inject_sources(model, si_cat, **kwargs):
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
        new_model = inject_sources_into_l2(model, si_cat, psftype="epsf", **kwargs)
    elif isinstance(model, MosaicModel):
        #  inject_sources_into_l3
        new_model = inject_sources_into_l3(model, si_cat, psftype="epsf", **kwargs)
    else:
        raise ValueError("The input model must be an ImageModel or MosaicModel.")

    return new_model


def make_cosmoslike_catalog(cen, ra, dec, exptimes, filters=None, seed=50, **kwargs):
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
    cen : astropy.coordinates.SkyCoord
        Location around which to generate sources.
    ra, dec : array_like (float)
        ra, dec of each source in objlist (degrees)
    exptimes : dictionary
        Exposure time for each filter (key)
    filters : list of str
        List of filter names
    seed : int
        Seed for random number generator

    Returns
    -------
    all_cat : astropy.Table
        Table for use with table_to_catalog to generate catalog for simulation.
    """
    from romanisim import bandpass, catalog

    # WFI bandpasses
    if filters is None:
        filters = set(bandpass.galsim2roman_bandpass.values())
    # Ensure the J-band is included, as it is used for sizes
    elif "F129" not in filters:
        filters.append("F129")

    # Set random source index for the catalog
    rng_numpy = np.random.default_rng(seed)

    # 75% of objects as galaxies, 25% of objects as point sources
    num_stars = int(len(ra) / 4)
    num_gals = len(ra) - num_stars

    # Galaxies
    # RomanISim's make cosmos galaxies method will return cosmos like objects with
    # uniformly distributed position angles. Radius and area are set to return
    # a full sample (~346k galaxies)
    gal_cat = catalog.make_cosmos_galaxies(
        cen, radius=1.0, bandpasses=filters, cat_area=(np.pi), **kwargs
    )

    # Trim to the required number of objects
    gal_cat = gal_cat[:num_gals]

    # Set magnitude spread
    mags = rng_numpy.uniform(low=-4.0, high=1.0, size=num_gals)

    # Brightest color flux
    gal_cat_data = (
        gal_cat.as_array(names=filters)
        .view(dtype=float)
        .reshape((len(gal_cat), len(filters)))
    )
    max_flux = gal_cat_data.max(axis=1)

    # Find gal mag limit
    gal_mag_limit = []
    for bp in filters:
        # Normalize the mag limit to exptimes
        if bp in exptimes:
            gal_mag_limit.append(
                HRGALMAGLIMIT[bp]
                + (1.25 * np.log10((exptimes[bp] * u.s).to(u.hour).value))
            )
    mag_tot = mags + max(gal_mag_limit)

    # Set bandpass fluxes
    for bp in filters:
        # Set scaled fluxes
        # flux_band = (flux_band / max_flux) * flux_spread_about_mag_limit
        gal_cat[bp] = ((gal_cat[bp] / max_flux) * (10.0 ** (-(mag_tot) / 2.5))).astype(
            "f4"
        )

    # Sizes are drawn from a log-normal distribution of J-band magnitudes
    # The parameters below are derived from the distribution for sizes
    # mu = (-0.1555 * (J - 17)) - 3.55
    jmag = -2.5 * np.log10(
        gal_cat["F129"],
        where=(gal_cat["F129"] > 0),
        out=np.array([-20.0] * len(gal_cat)),
    )
    mu = (-0.1555 * (jmag - 17)) - 3.55
    sigma = 0.15
    radsize = rng_numpy.normal(mu, sigma)
    gal_cat["half_light_radius"] = (10**radsize) * 3600 * u.arcsec

    # Set the radius floor
    gal_cat["half_light_radius"][gal_cat["half_light_radius"] < 0.036] = (
        0.036 * u.arcsec
    )

    # Randomize concentrations
    gal_cat["n"] = rng_numpy.uniform(low=0.8, high=4.5, size=num_gals)

    # Stars
    # Create base table
    star_cat = table.Table()
    star_cat["ra"] = num_stars * [0]
    star_cat["dec"] = num_stars * [0]
    star_cat["type"] = num_stars * ["PSF"]
    star_cat["n"] = num_stars * [-1]
    star_cat["half_light_radius"] = num_stars * [0] * u.arcsec
    star_cat["pa"] = num_stars * [0]
    star_cat["ba"] = num_stars * [1]

    # Set magnitude spread
    mags = rng_numpy.uniform(low=-4.0, high=1.0, size=num_stars)

    # Find point mag limit
    point_band_mag_limit = []
    for bp in filters:
        # Normalize the mag limit to exptimes
        if bp in exptimes:
            point_band_mag_limit.append(
                HRPOINTMAGLIMIT[bp]
                + (1.25 * np.log10((exptimes[bp] * u.s).to(u.hour).value))
            )
    point_mag_limit = max(point_band_mag_limit)

    # Set bandpass fluxes
    for bp in filters:
        star_cat[bp] = (10.0 ** (-(mags + point_mag_limit) / 2.5)).astype("f4")

    # Combine the objects
    all_cat = table.vstack([gal_cat, star_cat])

    # Set the positions to randomly selected objects
    all_cat["ra"] = np.array(ra).tolist()
    all_cat["dec"] = np.array(dec).tolist()

    return all_cat


def make_source_grid(
    model, yxmax=(5000, 5000), yxoffset=(50, 50), yxgrid=(20, 20), seed=50, **kwargs
):
    """
    Generate a grid of points to inject sources onto. The grid is set to the yxmax
    size for consistent spacing across input, and cropped as needed. An edge offset
    is specified within the parameters to avoid sources in that area.

    Parameters
    ----------
    model : `ImageModel` or `MosaicModel`
        Model into which to inject sources.
    yxmax : tuple of two ints
        Maximum extend of the grid
    yxoffset : int or tuple of two ints
        Edge Offset to place grid within
    yxgrid : tuple of two ints
        Grid point dimensions
    seed : int
        Seed for random number generator

    Returns
    -------
    y_pos_idx, x_pos_idx : array_like (int)
        y, x positions of each valid grid point yxmax
    """
    # Set random source index for the grid
    rng_numpy = np.random.default_rng(seed)

    if isinstance(yxoffset, (int, float)):
        yxoffset = (yxoffset, yxoffset)

    yspread, xspread = np.subtract((yxmax), 2 * np.array(yxoffset))
    yspace, xspace = np.ceil(np.divide((yspread, xspread), yxgrid))

    y0, x0 = (
        yxoffset[0] + rng_numpy.uniform(high=yspace),
        yxoffset[1] + rng_numpy.uniform(high=xspace),
    )

    # ypts, xpts = np.arange(yspace), np.arange(xspace)
    ypts, xpts = (
        np.arange(yxgrid[0], dtype=np.float64),
        np.arange(yxgrid[1], dtype=np.float64),
    )

    ypts *= yspace
    ypts += y0
    ypts += rng_numpy.uniform(low=-0.5, high=0.5, size=len(ypts))

    xpts *= xspace
    xpts += x0
    xpts += rng_numpy.uniform(low=-0.5, high=0.5, size=len(xpts))

    # Discard off-image positions
    ypts = ypts[ypts < (model.data.shape[0] - int(yxoffset[0]))]
    xpts = xpts[xpts < (model.data.shape[1] - int(yxoffset[1]))]

    # Create grid
    y_pos, x_pos = np.meshgrid(ypts, xpts)
    y_pos = np.ravel(y_pos)
    x_pos = np.ravel(x_pos)
    y_pos_idx, x_pos_idx = y_pos.astype(int), x_pos.astype(int)

    # Discard positions in NA empty regions
    nanmask = np.isnan(model.data[y_pos_idx, x_pos_idx])
    y_pos_idx, x_pos_idx = y_pos_idx[~nanmask], x_pos_idx[~nanmask]

    return y_pos_idx, x_pos_idx
