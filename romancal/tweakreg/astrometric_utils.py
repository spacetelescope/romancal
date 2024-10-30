import os

import requests
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord

from ..assign_wcs import utils as wcsutil

ASTROMETRIC_CAT_ENVVAR = "ASTROMETRIC_CATALOG_URL"
DEF_CAT_URL = "http://gsss.stsci.edu/webservices"

# VO request timeout (in seconds)
TIMEOUT = 30.0

if ASTROMETRIC_CAT_ENVVAR in os.environ:
    SERVICELOCATION = os.environ[ASTROMETRIC_CAT_ENVVAR]
else:
    SERVICELOCATION = DEF_CAT_URL


"""

Utility functions for creating an astrometric reference catalog.

"""


def compute_radius(wcs):
    """Compute the radius from the center to the furthest edge of the WCS."""

    fiducial = wcsutil.compute_fiducial([wcs], wcs.bounding_box)
    img_center = SkyCoord(ra=fiducial[0] * u.degree, dec=fiducial[1] * u.degree)
    wcs_foot = wcs.footprint()
    img_corners = SkyCoord(ra=wcs_foot[:, 0] * u.degree, dec=wcs_foot[:, 1] * u.degree)
    radius = img_center.separation(img_corners).max().value

    return radius, fiducial


def get_catalog(ra, dec, epoch=2016.0, sr=0.1, catalog="GAIADR3", timeout=TIMEOUT):
    """Extract catalog from VO web service.

    Parameters
    ----------
    ra : float
        Right Ascension (RA) of center of field-of-view (in decimal degrees)

    dec : float
        Declination (Dec) of center of field-of-view (in decimal degrees)

    epoch: float, optional
        Reference epoch used to update the coordinates for proper motion
        (in decimal year). Default: 2016.0.

    sr : float, optional
        Search radius (in decimal degrees) from field-of-view center to use
        for sources from catalog.  Default: 0.1 degrees

    catalog : str, optional
        Name of catalog to query, as defined by web-service.  Default: 'GAIADR3'

    timeout : float, optional
        Set the request timeout (in seconds). Default: 30 s.

    Returns
    -------
        A Table object of returned sources with all columns as provided by catalog.

    """
    service_type = "vo/CatalogSearch.aspx"
    spec_str = "RA={}&DEC={}&EPOCH={}&SR={}&FORMAT={}&CAT={}&MINDET=5"
    headers = {"Content-Type": "text/csv"}
    fmt = "CSV"

    spec = spec_str.format(ra, dec, epoch, sr, fmt, catalog)
    service_url = f"{SERVICELOCATION}/{service_type}?{spec}"
    try:
        rawcat = requests.get(service_url, headers=headers, timeout=timeout)
    except requests.exceptions.ConnectionError:
        raise requests.exceptions.ConnectionError(
            "Could not connect to the VO API server. Try again later."
        )
    except requests.exceptions.Timeout:
        raise requests.exceptions.Timeout("The request to the VO API server timed out.")
    except requests.exceptions.RequestException:
        raise requests.exceptions.RequestException(
            "There was an unexpected error with the request."
        )
    # convert from bytes to a String
    r_contents = rawcat.content.decode()
    rstr = r_contents.split("\r\n")
    # remove initial line describing the number of sources returned
    # CRITICAL to proper interpretation of CSV data
    del rstr[0]
    if len(rstr) == 0:
        raise Exception(
            """VO catalog service returned no results.\n
            Hint: maybe reviewing the search parameters might help."""
        )

    return table.Table.read(rstr, format="csv")
