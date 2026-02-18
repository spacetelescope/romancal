"""
This module determines which sky cells overlap with the given image.

Currently this assumes that the sky projected borders of all calibrated L2 images are straight.
"""

import logging
from functools import cached_property

import numpy as np
import spherical_geometry.great_circle_arc as sga
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from gwcs import WCS
from numpy.typing import NDArray

import romancal.skycell.skymap as sc

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ImageFootprint", "find_skycell_matches"]


class ImageFootprint:
    """abstraction of an image footprint"""

    _radec_corners: NDArray[float]

    def __init__(self, radec_vertices: list[tuple[float, float]]):
        """
        Parameters
        ----------
        radec_vertices: list[tuple[float, float]]
            vertices (usually the corners) of the image in right ascension and declination
        """
        self._radec_vertices = np.array(radec_vertices)

    @classmethod
    def from_wcs(
        cls,
        wcs: WCS,
        extra_vertices_per_edge: int = 0,
    ) -> "ImageFootprint":
        """create an image footprint from a GWCS object (and image shape, if no bounding box is present)

        Parameters
        ----------
        wcs: WCS :
            WCS object
        extra_vertices_per_edge: int :
            extra vertices to create on each edge to capture distortion (Default value = 0)

        Returns
        -------
        image footprint object
        """

        if (
            extra_vertices_per_edge <= 0
            and hasattr(wcs, "bounding_box")
            and wcs.bounding_box is not None
        ):
            vertex_points = wcs.footprint(center=False)
        else:
            vertex_points = np.array(
                sgp.SingleSphericalPolygon.from_wcs(
                    wcs, steps=extra_vertices_per_edge + 1
                ).to_lonlat()
            ).T

        return cls(vertex_points)

    @property
    def radec_corners(self) -> NDArray:
        """vertices in right ascension and declination in counterclockwise order"""
        return self._radec_vertices

    @cached_property
    def radec_center(self) -> tuple[float, float]:
        """center point in right ascension and declination"""
        return sgv.vector_to_lonlat(*self.vectorpoint_center)

    @cached_property
    def vectorpoint_vertices(self) -> NDArray[float]:
        """vertices in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(
            np.stack(sgv.lonlat_to_vector(*self.radec_corners.T), axis=1)
        )

    @cached_property
    def vectorpoint_center(self) -> tuple[float, float, float]:
        """center in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(np.mean(self.vectorpoint_vertices, axis=0))

    @cached_property
    def length(self) -> float:
        """diagonal length of the rectangular footprint"""
        # assume equally-spaced points around the perimeter
        # NOTE: this will produce an incorrect value with no error if the points are not equally spaced
        half_index_length = round(len(self.vectorpoint_vertices) / 2)
        return max(
            sga.length(
                self.vectorpoint_vertices[index],
                self.vectorpoint_vertices[index + half_index_length],
            )
            for index in range(len(self.vectorpoint_vertices) - half_index_length - 1)
        )

    @cached_property
    def circumference(self) -> float:
        """circumference of the rectangular footprint"""
        return sum(
            sga.length(
                self.vectorpoint_vertices[index], self.vectorpoint_vertices[index + 1]
            )
            for index in range(-1, len(self.vectorpoint_vertices) - 1)
        )

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this image footprint"""
        return sgp.SingleSphericalPolygon(
            points=self.vectorpoint_vertices,
            inside=self.vectorpoint_center,
        )

    @cached_property
    def area(self) -> float:
        """area of this footprint on the sphere in degrees squared"""
        return self.polygon.area()

    @cached_property
    def possibly_intersecting_projregions(self) -> int:
        """number of possibly intersecting projection regions"""
        if self.area > sc.SkyCells.area:
            return (
                # the number of times the smallest projection region could fit in the image footprint
                np.ceil(self.area / sc.ProjectionRegion.MIN_AREA)
                # plus multiplier for partial intersections on the perimeter
                * 8
            )
        else:
            # 4 foundational intersections
            return 4

    @cached_property
    def possibly_intersecting_skycells(self) -> int:
        """number of possibly intersecting skycells"""
        if self.polygon.area() > sc.SkyCells.area:
            return (
                # number of times a skycell could fit in the image footprint
                np.ceil(self.polygon.area() / sc.SkyCells.area)
                # plus multiplier for partial intersections on the perimeter
                * 8
            )
        else:
            # 4 foundational intersections
            return 4

    @cached_property
    def possible_intersecting_projregion_distance(self) -> float:
        """maximum possible distance to the center of an intersecting projection region"""
        return (self.length + sc.ProjectionRegion.MAX_LENGTH) / 2.0

    @cached_property
    def possible_intersecting_skycell_distance(self) -> float:
        """maximum possible distance to the center of an intersecting sky cell"""
        return (self.length + sc.SkyCells.length) / 2.0

    def __str__(self) -> str:
        return f"footprint {self.radec_corners}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.radec_corners!r})"


def find_skycell_matches(
    image_corners: list[tuple[float, float]] | NDArray[float] | WCS,
    skymap: sc.SkyMap = None,
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

    Returns
    -------
    Indices of all skycells (from the loaded skymap reference file) that overlap the supplied image.
    """

    if isinstance(image_corners, WCS):
        footprint = ImageFootprint.from_wcs(image_corners, extra_vertices_per_edge=3)
    else:
        footprint = ImageFootprint(image_corners)

    if skymap is None:
        skymap = sc.SKYMAP

    intersecting_skycell_indices = []

    # query the global k-d tree of projection regions for possible intersection candidates in (normalized) 3D space
    nearby_projregion_indices = np.array(
        skymap.projection_regions_kdtree.query_ball_point(
            footprint.vectorpoint_center,
            r=footprint.possible_intersecting_projregion_distance * 1.1,
        )
    )
    nearby_projregion_indices = nearby_projregion_indices[
        nearby_projregion_indices != len(skymap.model.projection_regions)
    ]

    for projregion_index in nearby_projregion_indices:
        projregion = sc.ProjectionRegion(projregion_index)
        if footprint.polygon.intersects_poly(projregion.polygon):
            # query the LOCAL k-d tree of skycells for possible intersection candidates in (normalized) 3D space
            projregion_nearby_skycell_indices = np.array(
                projregion.skycells.kdtree.query_ball_point(
                    footprint.vectorpoint_center,
                    r=footprint.possible_intersecting_skycell_distance * 1.1,
                )
            )

            projregion_nearby_skycells = sc.SkyCells(
                np.array(
                    projregion_nearby_skycell_indices[
                        projregion_nearby_skycell_indices != len(projregion.skycells)
                    ]
                )
                + projregion.data["skycell_start"],
                skymap=skymap,
            )

            # find polygons that intersect the image footprint
            for skycell_index, skycell_polygon in zip(
                projregion_nearby_skycells.indices,
                projregion_nearby_skycells.polygons,
                strict=True,
            ):
                if footprint.polygon.intersects_poly(skycell_polygon):
                    intersecting_skycell_indices.append(skycell_index)

    return intersecting_skycell_indices
