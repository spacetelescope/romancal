"""Test utilities"""

import numpy as np
from astropy.stats import mad_std


def comp_wcs_grids_arcs(wcs_a, wcs_b, npix=4088, interval=10):
    """Compare world grids produced by the two wcs

    Parameters
    ----------
    wcs_a, wcs_b : gwcs.WCS
        The wcs object to compare.

    npix : int
        The size of the grid to produce.

    interval : int
        The interval to check over.

    Returns
    -------
    mad_std : float
        The numpy mad_std in mas
    """
    xx, yy = np.meshgrid(np.linspace(0, npix, interval), np.linspace(0, npix, interval))
    ra_a, dec_a = wcs_a(xx, yy, with_bounding_box=False)
    ra_b, dec_b = wcs_b(xx, yy, with_bounding_box=False)
    dec_med = np.nanmedian(dec_b)

    ra_mad = mad_std((ra_a - ra_b) * np.cos(np.radians(dec_med)), ignore_nan=True
                     ) * 60.0 * 60.0 * 1000.0
    dec_mad = mad_std(dec_a - dec_b, ignore_nan=True) * 60.0 * 60.0 * 1000.0

    return ra_mad, dec_mad
