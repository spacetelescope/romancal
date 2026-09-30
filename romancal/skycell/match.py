"""
This module determines which sky cells overlap with the given image.

Matching proceeds in three steps:

1. keep projection regions whose skycells could reach the image footprint,
   judged by the angular distance between their centers;
2. within those regions, keep skycells whose centers could lie close enough
   to the footprint to overlap it;
3. test the remaining skycells for overlap exactly, in the tangent plane of
   their projection region.

Each skycell is grown by a small buffer before the exact test, so that a
skycell is kept if it comes within the buffer of the footprint.  This makes
matching conservative: at worst a few extra skycells are selected.  A
non-convex footprint is replaced by its convex hull, which is conservative in
the same way.

This assumes that the sky projected borders of all calibrated L2
images are great circles.  The extra buffer handles the deviation from
great circles.
"""

import logging
from functools import cached_property

import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from astropy.modeling import Model, models
from gwcs import WCS
from numpy.typing import NDArray
from scipy.spatial import ConvexHull

import romancal.skycell.skymap as sc

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["find_skycell_matches"]


class _ImageFootprint:
    """abstraction of an image footprint"""

    _radec_vertices: NDArray[float]

    def __init__(self, radec_vertices: list[tuple[float, float]]):
        """
        Parameters
        ----------
        radec_vertices: list[tuple[float, float]]
            vertices (usually the corners) of the image in right ascension and declination
        """
        self._radec_vertices = np.array(radec_vertices)

    @classmethod
    def from_wcs(cls, wcs: WCS) -> "_ImageFootprint":
        """create an image footprint from the corners of a GWCS object

        Also uses image shape, if no bounding box is present

        Parameters
        ----------
        wcs: WCS :
            WCS object

        Returns
        -------
        image footprint object
        """

        if hasattr(wcs, "bounding_box") and wcs.bounding_box is not None:
            vertex_points = wcs.footprint(center=False)
        else:
            # the polygon is closed, repeating its first vertex at the end
            vertex_points = np.array(
                sgp.SingleSphericalPolygon.from_wcs(wcs, steps=1).to_lonlat()
            ).T[:-1]

        return cls(vertex_points)

    @property
    def radec_corners(self) -> NDArray:
        """vertices in right ascension and declination"""
        return self._radec_vertices

    @cached_property
    def vectorpoint_vertices(self) -> NDArray[float]:
        """vertices in 3D Cartesian space on the unit sphere"""
        return sc._vectorpoints(self.radec_corners)

    @cached_property
    def vectorpoint_center(self) -> NDArray[float]:
        """center in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(np.mean(self.vectorpoint_vertices, axis=0))

    @cached_property
    def radius(self) -> float:
        """largest angular distance in degrees from the center to any point in the footprint"""
        # the farthest point of a small spherical polygon is one of its vertices
        return sc._separation(self.vectorpoint_vertices, self.vectorpoint_center).max()

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this image footprint"""
        return sgp.SingleSphericalPolygon(
            points=self.vectorpoint_vertices,
            inside=self.vectorpoint_center,
        )

    def __str__(self) -> str:
        return f"footprint {self.radec_corners}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.radec_corners!r})"


def find_skycell_matches(
    image_corners: list[tuple[float, float]] | NDArray[float] | WCS,
    skymap: sc.SkyMap = None,
    buffer_pixels: float = 20,
) -> list[int]:
    """Find sky cells overlapping the provided image footprint

    Parameters
    ----------
    image_corners : list | np.ndarray | WCS :
        Either a squence of 4 (ra, dec) pairs, or
        equivalent 2-d numpy array, or a GWCS instance.
        A GWCS instance must have `.bounding_box` or `.pixel_shape` attribute defined.
    skymap : sc.SkyMap :
        skymap instance; defaults to global SKYMAP (Default value = None)
    buffer_pixels : float :
        also match sky cells that come within this many skycell pixels of the
        image footprint, to allow for distortion of the image edges
        (Default value = 20)

    Returns
    -------
    Indices of all skycells (from the loaded skymap reference file) that overlap the supplied image.
    """

    if isinstance(image_corners, WCS):
        footprint = _ImageFootprint.from_wcs(image_corners)
    else:
        footprint = _ImageFootprint(image_corners)

    if skymap is None:
        skymap = sc.SKYMAP

    # all angles in degrees
    buffer = buffer_pixels * skymap.pixel_scale
    # tangent plane projection shrinks angles, so no point in a skycell is
    # farther from its center than half the diagonal of its pixel grid
    skycell_radius = skymap.pixel_shape[0] / np.sqrt(2) * skymap.pixel_scale

    # 1. projection regions whose skycells could reach the footprint
    projregion_separations = sc._separation(
        skymap._projection_region_vectorpoints, footprint.vectorpoint_center
    )
    nearby_projregion_indices = np.flatnonzero(
        projregion_separations
        <= skymap._projection_region_radii + footprint.radius + buffer
    )

    intersecting_skycell_indices = []
    for projregion_index in nearby_projregion_indices:
        projregion = sc.ProjectionRegion(projregion_index, skymap=skymap)

        # 2. skycells in this region whose centers are near the footprint
        skycells = sc.SkyCells(projregion.skycell_indices, skymap=skymap)
        skycell_separations = sc._separation(
            sc._vectorpoints(skycells.radec_centers), footprint.vectorpoint_center
        )
        nearby = skycell_separations <= footprint.radius + skycell_radius + buffer
        skycells = sc.SkyCells(skycells.indices[nearby], skymap=skymap)

        # 3. skycells that intersect the footprint
        intersecting_skycell_indices.extend(
            skycells.indices[
                _intersects(footprint, skycells, projregion, buffer)
            ].tolist()
        )

    return intersecting_skycell_indices


def _tangent_plane(ra: float, dec: float) -> Model:
    """gnomonic (tangent plane) projection at the given point, in degrees"""
    return models.RotateCelestial2Native(ra, dec, 180) | models.Sky2Pix_TAN()


def _intersects(
    footprint: _ImageFootprint,
    skycells: sc.SkyCells,
    projregion: sc.ProjectionRegion,
    buffer: float,
) -> NDArray[bool]:
    """whether each skycell, grown by `buffer` degrees, intersects the footprint

    This routine takes advantage of the fact that the skycell tessellation
    has straight line boundaries in the tangent plane projections of each
    projection region, and that the great circle edges of the image footprint
    correspond to straight lines in the tangent plane projection.

    Because of this we can test for overlap using 2D polygons without
    spherical geometry.
    """
    projection = _tangent_plane(*projregion.radec_tangent)

    # footprint vertices in the tangent plane: (vertex, x or y)
    polygon = np.stack(projection(*footprint.radec_corners.T), axis=-1)
    if not np.all(np.isfinite(polygon)):
        raise ValueError("footprint is too large to project onto a tangent plane")

    # take the convex hull to simplify overlap logic; practically most
    # footprints will be convex anyway.  Its vertices are in order around it.
    polygon = polygon[ConvexHull(polygon).vertices]

    # skycell corners in the tangent plane: (skycell, corner, x or y)
    rectangles = np.stack(
        projection(skycells.radec_corners[..., 0], skycells.radec_corners[..., 1]),
        axis=-1,
    )
    # grow each rectangle by moving its corners away from its center
    centers = rectangles.mean(axis=1, keepdims=True)
    rectangles = rectangles + buffer * np.sign(rectangles - centers)
    # (skycell, x or y): the (left, bottom) and (right, top) of each rectangle
    lower, upper = rectangles.min(axis=1), rectangles.max(axis=1)

    # A rectangle and a convex polygon are disjoint if and only if the line
    # along one of their edges separates them.  First, the rectangle edges:
    # is the whole polygon left of the rectangle's left edge, or below its
    # bottom edge, or right of its right edge, or above its top edge?
    polygon_lower, polygon_upper = polygon.min(axis=0), polygon.max(axis=0)
    beyond_edge = (polygon_upper < lower) | (polygon_lower > upper)
    separated = np.any(beyond_edge, axis=1)

    # then the polygon edges: do all rectangle corners lie strictly on the
    # other side of the edge from the polygon?
    def side(start, end, points):
        """which side of the line from `start` to `end` each point lies on"""
        edge, offset = end - start, points - start
        return np.sign(edge[0] * offset[..., 1] - edge[1] * offset[..., 0])

    for start, end in zip(polygon, np.roll(polygon, -1, axis=0), strict=True):
        polygon_side = side(start, end, polygon.mean(axis=0))
        separated |= np.all(side(start, end, rectangles) == -polygon_side, axis=1)

    return ~separated
