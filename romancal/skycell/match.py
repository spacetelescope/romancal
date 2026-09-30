"""
This module determines which sky cells overlap with the given image.

Matching proceeds in three steps:

1. keep projection regions whose skycells could reach the image footprint,
   judged by the angular distance between their centers;
2. within those regions, keep skycells whose centers could lie close enough
   to the footprint to overlap it;
3. test the remaining skycells for overlap exactly.

Each skycell is grown by a small buffer before the exact test, so that a
skycell is kept if it comes within the buffer of the footprint.  This makes
matching conservative: at worst a few extra skycells are selected.

Currently this assumes that the sky projected borders of all calibrated L2
images are great circles; the buffer is meant to cover the difference.
"""

import logging
from functools import cached_property

import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from gwcs import WCS
from numpy.typing import NDArray

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
        """create an image footprint from the corners of a GWCS object (and image shape, if no bounding box is present)

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
        return np.stack(sgv.lonlat_to_vector(*self.radec_corners.T), axis=1)

    @cached_property
    def vectorpoint_center(self) -> NDArray[float]:
        """center in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(np.mean(self.vectorpoint_vertices, axis=0))

    @cached_property
    def radius(self) -> float:
        """largest angular distance in radians from the center to any point in the footprint"""
        # the farthest point of a small spherical polygon is one of its vertices
        return _separation(self.vectorpoint_vertices, self.vectorpoint_center).max()

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

    pixel_scale = np.radians(skymap.pixel_scale)
    buffer = buffer_pixels * pixel_scale
    # gnomonic projection does not stretch angles, so no point in a skycell is
    # farther from its center than half the diagonal of its pixel grid
    skycell_radius = skymap.pixel_shape[0] / np.sqrt(2) * pixel_scale

    # 1. projection regions whose skycells could reach the footprint
    nearby_projregion_indices = np.nonzero(
        _separation(
            skymap._projection_region_vectorpoints, footprint.vectorpoint_center
        )
        <= skymap._projection_region_radii + footprint.radius + buffer
    )[0]

    intersecting_skycell_indices = []
    for projregion_index in nearby_projregion_indices:
        projregion = sc.ProjectionRegion(projregion_index, skymap=skymap)

        # 2. skycells in this region whose centers are near the footprint
        skycells = sc.SkyCells(projregion.skycell_indices, skymap=skymap)
        nearby = (
            _separation(
                _vectorpoints(skycells.radec_centers), footprint.vectorpoint_center
            )
            <= footprint.radius + skycell_radius + buffer
        )
        skycells = sc.SkyCells(skycells.indices[nearby], skymap=skymap)

        # 3. skycells that intersect the footprint
        intersecting_skycell_indices.extend(
            skycells.indices[_intersects(footprint, skycells, buffer)].tolist()
        )

    return intersecting_skycell_indices


def _vectorpoints(radec: NDArray[float]) -> NDArray[float]:
    """convert (..., 2) right ascension and declination to (..., 3) unit vectors"""
    return np.stack(sgv.lonlat_to_vector(radec[..., 0], radec[..., 1]), axis=-1)


def _separation(vectorpoints: NDArray[float], vectorpoint: NDArray[float]):
    """angular distance in radians between unit vectors"""
    return np.arccos(np.clip(np.sum(vectorpoints * vectorpoint, axis=-1), -1, 1))


def _intersects(
    footprint: _ImageFootprint, skycells: sc.SkyCells, buffer: float
) -> NDArray[bool]:
    """whether each skycell, grown by `buffer` radians, intersects the footprint"""
    centers = _vectorpoints(skycells.radec_centers)[:, None, :]
    corners = _vectorpoints(skycells.radec_corners)
    if buffer > 0:
        # move each corner away from the skycell center; moving the corners of
        # a square by buffer * sqrt(2) moves its edges out by buffer
        cos_angle = np.clip(np.sum(corners * centers, axis=-1, keepdims=True), -1, 1)
        direction = corners - cos_angle * centers
        direction /= np.linalg.norm(direction, axis=-1, keepdims=True)
        angle = np.arccos(cos_angle) + buffer * np.sqrt(2)
        corners = np.cos(angle) * centers + np.sin(angle) * direction

    return np.array(
        [
            footprint.polygon.intersects_poly(
                sgp.SingleSphericalPolygon(points=points, inside=center)
            )
            for points, center in zip(corners, centers[:, 0], strict=True)
        ],
        dtype=bool,
    )
