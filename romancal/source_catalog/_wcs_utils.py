"""
Local WCS geometry: the pixel solid angle and the direction of celestial
North, sampled across an image.
"""

import astropy.units as u
import numpy as np
from scipy.interpolate import RectBivariateSpline


def _wrap_lon_diff(lon1, lon0):
    """
    Difference two longitudes (in degrees), wrapping across the 0/360
    degree branch cut.
    """
    return (lon1 - lon0 + 180.0) % 360.0 - 180.0


def wcs_jacobian(wcs, x, y):
    """
    Compute the local WCS Jacobian at the given pixel positions.

    The Jacobian is evaluated by central differences over one pixel and
    describes the sky offset produced by a unit step along each detector
    axis::

        [[d(alpha cos delta)/dx, d(alpha cos delta)/dy],
         [d(delta)/dx,           d(delta)/dy]]

    where ``alpha`` and ``delta`` are the longitude and latitude of the
    world frame. The longitude derivatives carry the ``cos(delta)``
    factor so that both rows are true angles on the sky, in arcsec.
    Everything about the local pixel geometry follows from this matrix:
    its determinant is the pixel solid angle, and the direction mapping
    to increasing ``delta`` is celestial North.

    Parameters
    ----------
    wcs : WCS object
        A world coordinate system transformation mapping pixel to world
        (longitude, latitude) coordinates in degrees.

    x, y : `~numpy.ndarray`
        Pixel coordinates, of any common shape.

    Returns
    -------
    jacobian : `~numpy.ndarray`
        Array of shape ``x.shape + (2, 2)`` in arcsec per pixel.
    """

    # Positions may lie outside the WCS bounding box, which would
    # otherwise be evaluated as NaN.
    def evaluate(xx, yy):
        return wcs(xx, yy, with_bounding_box=False)

    _, lat = evaluate(x, y)
    cos_lat = np.cos(np.radians(lat))

    lon_xlo, lat_xlo = evaluate(x - 0.5, y)
    lon_xhi, lat_xhi = evaluate(x + 0.5, y)
    lon_ylo, lat_ylo = evaluate(x, y - 0.5)
    lon_yhi, lat_yhi = evaluate(x, y + 0.5)

    arcsec = 3600.0
    dxi_dx = _wrap_lon_diff(lon_xhi, lon_xlo) * cos_lat * arcsec
    dxi_dy = _wrap_lon_diff(lon_yhi, lon_ylo) * cos_lat * arcsec
    deta_dx = (lat_xhi - lat_xlo) * arcsec
    deta_dy = (lat_yhi - lat_ylo) * arcsec

    return np.stack(
        [np.stack([dxi_dx, dxi_dy], axis=-1), np.stack([deta_dx, deta_dy], axis=-1)],
        axis=-2,
    )


def pixel_area_map(wcs, shape, step=64):
    """
    Compute the on-sky solid angle of every pixel in an image.

    The area is the determinant of the local WCS Jacobian, evaluated on
    a coarse grid and then spline-interpolated onto the full image grid.
    The Roman distortion varies smoothly on scales far larger than
    ``step``, so the interpolation error (~3e-5 at the default ``step``
    for a WFI detector) is negligible compared to the ~2.5% peak-to-peak
    area variation the map exists to capture.

    Parameters
    ----------
    wcs : WCS object
        A world coordinate system transformation mapping pixel to world
        (longitude, latitude) coordinates in degrees.

    shape : tuple of int
        The ``(ny, nx)`` shape of the image.

    step : int, optional
        Spacing, in pixels, of the coarse grid on which the Jacobian is
        evaluated.

    Returns
    -------
    area : `~astropy.units.Quantity`
        Array of shape ``shape`` giving the solid angle of each pixel in
        steradians.
    """
    ny, nx = shape

    # Pad the coarse grid by one step so that the interpolation covers
    # the full image without extrapolating.
    gy = np.arange(-step, ny + step, step, dtype=float)
    gx = np.arange(-step, nx + step, step, dtype=float)
    xx, yy = np.meshgrid(gx, gy)

    coarse_area = np.abs(np.linalg.det(wcs_jacobian(wcs, xx, yy)))

    spline = RectBivariateSpline(gy, gx, coarse_area)
    area = spline(np.arange(ny, dtype=float), np.arange(nx, dtype=float))

    return (area.astype(np.float32) * u.arcsec**2).to(u.sr)


def pixel_area_at(area_map, x, y):
    """
    Look up the pixel solid angle at each of the given positions.

    Parameters
    ----------
    area_map : `~astropy.units.Quantity`
        2D array of per-pixel solid angles, as returned by
        `pixel_area_map`.

    x, y : `~numpy.ndarray`
        Pixel coordinates. Non-finite positions fall back to the first
        pixel; the other properties of such sources are non-finite in
        any case.

    Returns
    -------
    area : `~astropy.units.Quantity`
        The solid angle of the pixel containing each position.
    """
    ny, nx = area_map.shape
    x = np.clip(np.round(np.nan_to_num(x)).astype(int), 0, nx - 1)
    y = np.clip(np.round(np.nan_to_num(y)).astype(int), 0, ny - 1)
    return area_map[y, x]


def north_angle_at(wcs, x, y):
    """
    The position angle of celestial North at each of the given
    positions.

    The angle is measured counterclockwise from the positive x axis of
    the detector, matching the convention of the image-center value it
    replaces.

    North is the pixel-space direction that maps to a pure
    increase in latitude, so it is ``J^-1 [0, 1]``, using the Jacobian
    J.  This needs only the forward transform, unlike offsetting north
    on the sky and inverting, which fails near an array edge and at the
    pole.

    Only the direction of ``J^-1 [0, 1]`` matters here, and for a 2x2
    matrix that is ``[-J[0, 1], J[0, 0]]`` up to a factor of the
    determinant, so the determinant enters only through its sign.
    Avoiding the division keeps this finite at a celestial pole, where
    J is singular; North is genuinely undefined there, but returning an
    arbitrary angle is better than failing for every source in the
    image.

    Parameters
    ----------
    wcs : WCS object
        A world coordinate system transformation mapping pixel to world
        (longitude, latitude) coordinates in degrees.

    x, y : `~numpy.ndarray`
        Pixel coordinates, of any common shape.

    Returns
    -------
    angle : `~astropy.units.Quantity`
        The position angle of North at each position, in degrees.
    """
    jacobian = wcs_jacobian(wcs, np.asarray(x), np.asarray(y))
    sign = np.sign(np.linalg.det(jacobian))
    north_x = -sign * jacobian[..., 0, 1]
    north_y = sign * jacobian[..., 0, 0]
    return (np.degrees(np.arctan2(north_y, north_x)) * u.deg).astype(np.float32)
