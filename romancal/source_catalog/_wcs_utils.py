"""
Local WCS geometry. The pixel solid angle is sampled across an image.
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
    norm of their cross product is the pixel's solid angle.  We use 3D
    cartesian vectors rather than e.g. angular sky coordinates to avoid
    any singularities.

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


def pixel_area_from_wcs(wcs, x, y):
    """
    Compute the on-sky solid angle of the pixels at the given positions.

    The area is that of the parallelogram spanned by the pixel's two
    edge vectors (see `wcs_jacobian`).

    Parameters
    ----------
    wcs : WCS object
        A world coordinate system transformation mapping pixel to world
        (longitude, latitude) coordinates in degrees.

    x, y : `~numpy.ndarray`
        Pixel coordinates, of any common shape.

    Returns
    -------
    area : `~astropy.units.Quantity`
        Array of shape ``x.shape`` giving the solid angle of each pixel
        in arcsec**2.
    """
    jacobian = wcs_jacobian(wcs, x, y)
    area = np.linalg.norm(np.cross(jacobian[..., 0], jacobian[..., 1]), axis=-1)
    return area * u.arcsec**2


def pixel_to_sky_errors(wcs, x, y, x_err, y_err):
    """
    Convert independent pixel position uncertainties to uncertainties in sky coordinates.

    The position uncertainties in pixels are transformed to sky coordinates using the
    local WCS Jacobian (see `wcs_jacobian`).  Covariance between the x and
    y uncertainties is neglected.

    Parameters
    ----------
    wcs : WCS object
        A world coordinate system transformation mapping pixel to world
        coordinates in degrees.

    x, y : `~numpy.ndarray`
        Pixel coordinates.

    x_err, y_err : `~numpy.ndarray`
        The 1-sigma pixel position errors, with the same shape as ``x``.

    Returns
    -------
    east_err, north_err : `~astropy.units.Quantity`
        The 1-sigma errors in arcsec in the right ascension
        direction and declination directions.
    """
    # We want, e.g. for the East error,
    #
    #     east_err = sqrt((dE/dx * x_err)**2 + (dE/dy * y_err)**2)
    #
    # where dE/dx is the distance (arcsec) the sky position moves toward
    # East per pixel step along x. We get dE/dx in two steps: how the
    # position moves in 3D per pixel step (the Jacobian), and how much of
    # that motion points East (a dot product with the East unit vector).
    #
    # Shapes below are for N sources.

    # Step 1: the Jacobian, shape (N, 3, 2). jacobian[n, :, 0] is the
    # 3D (cartesian) step, in arcsec, of the sky position's unit vector
    # when moving one pixel along x; jacobian[n, :, 1] is the same for y.
    # Both steps lie in the plane tangent to the sky at the source.
    jacobian = wcs_jacobian(wcs, x, y)

    # Step 2: the East and North unit vectors in that tangent plane, as
    # 3D cartesian vectors of shape (N, 3). They are the derivatives of
    # the unit vector (cos(lat) cos(lon), cos(lat) sin(lon), sin(lat))
    # with respect to lon (divided by cos(lat)) and lat respectively.
    lon, lat = np.deg2rad(wcs(x, y, with_bounding_box=False))
    east = np.stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)], axis=-1)
    north = np.stack(
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)],
        axis=-1,
    )

    # Dot each pixel step with East and North, shape (N, 2):
    # d_east[n] = (dE/dx, dE/dy) and d_north[n] = (dN/dx, dN/dy), in
    # arcsec per pixel. Because East is a great-circle direction, dE
    # already includes the cos(dec) factor.
    d_east = np.einsum("ni,nij->nj", east, jacobian)
    d_north = np.einsum("ni,nij->nj", north, jacobian)

    east_err = np.hypot(d_east[..., 0] * x_err, d_east[..., 1] * y_err)
    north_err = np.hypot(d_north[..., 0] * x_err, d_north[..., 1] * y_err)
    return east_err * u.arcsec, north_err * u.arcsec


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

    coarse_area = pixel_area_from_wcs(wcs, xx, yy).to_value(u.arcsec**2)

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
