"""
Local WCS geometry: the pixel solid angle and the direction of celestial
North, sampled across an image.
"""

import astropy.units as u
import numpy as np
from astropy.coordinates import UnitSphericalRepresentation
from scipy.interpolate import RectBivariateSpline

ARCSEC_PER_RADIAN = (1.0 * u.rad).to_value(u.arcsec)


def _spherical(lon, lat):
    """
    A `~astropy.coordinates.UnitSphericalRepresentation` of the given
    sky positions, with longitude and latitude given in degrees.
    """
    return UnitSphericalRepresentation(lon * u.deg, lat * u.deg)


def _as_vectors(representation):
    """
    Reshape a `~astropy.coordinates.CartesianRepresentation` to put the
    Cartesian components last, giving shape ``lon.shape + (3,)``.
    """
    return np.moveaxis(representation.xyz.value, 0, -1)


def _direction_cosines(lon, lat):
    """
    The Cartesian unit vectors of the given sky positions.

    This representation is smooth over the whole sphere, unlike
    (longitude, latitude), which is degenerate at the poles and
    discontinuous at the 0/360 degree branch cut.
    """
    return _as_vectors(_spherical(lon, lat).to_cartesian())


def _north_vectors(lon, lat):
    """
    Unit vectors pointing toward celestial North at the given sky
    positions.

    This is the direction of increasing latitude, which
    `~astropy.coordinates.UnitSphericalRepresentation.unit_vectors`
    supplies directly. North is undefined at a pole; astropy resolves
    that by returning the limit reached along the given longitude,
    which is as good an answer as any and keeps the result finite.
    """
    return _as_vectors(_spherical(lon, lat).unit_vectors()["lat"])


def wcs_jacobian(wcs, x, y):
    """
    Compute the local WCS Jacobian at the given pixel positions.

    The Jacobian is evaluated by central differences over one pixel and
    describes the sky offset produced by a unit step along each detector
    axis. Its two columns are the pixel's edge vectors, expressed in the
    Cartesian direction cosines ``n`` of the sky position::

        [[dn_x/dx, dn_x/dy],
         [dn_y/dx, dn_y/dy],
         [dn_z/dx, dn_z/dy]]

    Differencing unit vectors rather than (longitude, latitude) keeps
    this well behaved everywhere: there is no ``cos(latitude)`` factor
    to collapse at a celestial pole, and no branch cut to wrap across at
    a longitude of zero. Everything about the local pixel geometry
    follows from these two vectors: the norm of their cross product is
    the pixel solid angle, and the pixel-space direction they map to
    celestial North gives the local position angle.

    The columns are chord lengths on the unit sphere rather than arc
    lengths, which understates the angle by a relative ``theta**2 / 24``:
    around 1e-14 for a pixel, and negligible beside the truncation error
    of the central difference itself.

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
        return _direction_cosines(*wcs(xx, yy, with_bounding_box=False))

    d_dx = (evaluate(x + 0.5, y) - evaluate(x - 0.5, y)) * ARCSEC_PER_RADIAN
    d_dy = (evaluate(x, y + 0.5) - evaluate(x, y - 0.5)) * ARCSEC_PER_RADIAN

    return np.stack([d_dx, d_dy], axis=-1)


def pixel_area_map(wcs, shape, step=64):
    """
    Compute the on-sky solid angle of every pixel in an image.

    The area is the area of the parallelogram spanned by the columns of
    the local WCS Jacobian, evaluated on a coarse grid and then
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
    replaces.

    The answer is the pixel-space step that the Jacobian J maps to the
    North direction, so it solves ``J v = north``. That system is
    overdetermined but consistent, since North lies in the plane the
    columns of J span, and it is solved in the least-squares sense
    through the 2x2 matrix ``J.T J``. Using only the forward transform
    this way avoids the array-edge and pole failures of offsetting north
    on the sky and inverting.

    Only the direction of ``v`` matters here, so the determinant of
    ``J.T J`` can be dropped: it is a Gram determinant and so never
    negative, and dividing by it would only rescale ``v``.

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
    x = np.asarray(x)
    y = np.asarray(y)
    jacobian = wcs_jacobian(wcs, x, y)
    north = _north_vectors(*wcs(x, y, with_bounding_box=False))

    gram = np.einsum("...ki,...kj->...ij", jacobian, jacobian)
    projection = np.einsum("...ki,...k->...i", jacobian, north)

    # The adjugate of the 2x2 Gram matrix, applied to the projection
    north_x = (
        gram[..., 1, 1] * projection[..., 0] - gram[..., 0, 1] * projection[..., 1]
    )
    north_y = (
        gram[..., 0, 0] * projection[..., 1] - gram[..., 1, 0] * projection[..., 0]
    )

    return (np.degrees(np.arctan2(north_y, north_x)) * u.deg).astype(np.float32)
