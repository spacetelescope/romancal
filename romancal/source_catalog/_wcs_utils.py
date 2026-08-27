"""
Local WCS geometry: the pixel solid angle and the direction of celestial
North, sampled across an image.
"""

import astropy.units as u
import numpy as np
from astropy.coordinates import UnitSphericalRepresentation
from scipy.interpolate import RectBivariateSpline

ARCSEC_PER_RADIAN = (1.0 * u.rad).to_value(u.arcsec)


def _unit_vectors(lon, lat):
    """
    The unit vectors pointing at the given sky positions, with shape
    ``lon.shape + (3,)``. Longitude and latitude are in degrees.

    Unit vectors are smooth over the whole sphere, unlike
    (longitude, latitude), which is degenerate at the poles and
    discontinuous at the 0/360 degree branch cut.
    """
    spherical = UnitSphericalRepresentation(lon * u.deg, lat * u.deg)
    return np.moveaxis(spherical.to_cartesian().xyz.value, 0, -1)


def wcs_jacobian(wcs, x, y):
    """
    Compute the local WCS Jacobian at the given pixel positions.

    The two columns of the Jacobian are the pixel's edge vectors: the
    step taken in xyz by moving one pixel along the detector x axis, and
    the same for the y axis. They are evaluated by central differences
    of the unit vector pointing at the sky position::

        [d(unit vector)/dx, d(unit vector)/dy]

    Everything about the local pixel geometry follows from these two
    vectors. They span the plane tangent to the sky at that pixel; the
    norm of their cross product is the pixel's solid angle, and the
    combination of them that reaches the celestial pole gives the local
    direction of North.  We use 3D cartesian vectors rather than
    e.g. angular sky coordinates to avoid any singularities.

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
        Array of shape ``x.shape + (3, 2)`` in arcsec per pixel.
    """

    # Positions may lie outside the WCS bounding box, which would
    # otherwise be evaluated as NaN.
    def evaluate(xx, yy):
        return _unit_vectors(*wcs(xx, yy, with_bounding_box=False))

    d_dx = (evaluate(x + 0.5, y) - evaluate(x - 0.5, y)) * ARCSEC_PER_RADIAN
    d_dy = (evaluate(x, y + 0.5) - evaluate(x, y - 0.5)) * ARCSEC_PER_RADIAN

    return np.stack([d_dx, d_dy], axis=-1)


def pixel_area_map(wcs, shape, step=64):
    """
    Compute the on-sky solid angle of every pixel in an image.

    The area is that of the parallelogram spanned by the pixel's two
    edge vectors, evaluated on a coarse grid and then
    spline-interpolated onto the full image grid. The Roman distortion
    varies smoothly on scales far larger than ``step``, so the
    interpolation error (~1e-7 at the default ``step`` for a WFI
    detector) is negligible compared to the ~2.5% peak-to-peak area
    variation the map exists to capture.

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

    jacobian = wcs_jacobian(wcs, xx, yy)
    coarse_area = np.linalg.norm(np.cross(jacobian[..., 0], jacobian[..., 1]), axis=-1)

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
    replaces and the convention `orientation_sky` expects.

    We use unit vectors on the sphere here to avoid degeneracies
    at the celestial pole.

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
    # pixel x and y edges in cartesian coordinates
    jacobian = wcs_jacobian(wcs, np.asarray(x), np.asarray(y))

    pole = np.array([0.0, 0.0, 1.0])[:, np.newaxis]
    transpose = jacobian.swapaxes(-1, -2)

    # step gives the combination of the edge vectors so that
    # pole ~ a * edge_x + b * edge_y
    # this is the least squares solution of J step ~ pole,
    # given by (JT J) step = JT pole
    # note: near the pole, most of the pole vector is out of the
    # plane of the pixel; JT limits to the projection along
    # the pixel edges so that we're solving only for the in-plane
    # component
    step = np.linalg.solve(transpose @ jacobian, transpose @ pole)
    north_x, north_y = step[..., 0, 0], step[..., 1, 0]

    return (np.degrees(np.arctan2(north_y, north_x)) * u.deg).astype(np.float32)
