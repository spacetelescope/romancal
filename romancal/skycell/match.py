"""
This module determines which sky cells overlap with the given image.

Currently this assumes that the sky projected borders of all calibrated L2 images are straight.
"""

import logging

import numpy as np
import sphersgeo
from gwcs import WCS
from numpy.typing import NDArray
from sphersgeo import ArcString, MultiSphericalPoint, SphericalPoint, SphericalPolygon

import romancal.skycell.skymap as sc

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ImageFootprint", "find_skycell_matches"]


class ImageFootprint:
    """abstraction of an image footprint"""

    def __init__(self, polygon: SphericalPolygon):
        self.polygon = polygon

    @classmethod
    def from_lonlats(
        cls, radec_vertices: list[tuple[float, float]]
    ) -> "ImageFootprint":
        """
        Parameters
        ----------
        radec_vertices: list[tuple[float, float]]
            vertices (usually the corners) of the image in right ascension and declination
        """
        return cls(
            SphericalPolygon(
                ArcString(MultiSphericalPoint.from_lonlats(radec_vertices))
            )
        )

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
        return cls(sphersgeo.from_wcs.polygon_from_wcs(wcs, extra_vertices_per_edge))

    @property
    def center(self) -> SphericalPoint:
        """center point of image footprint on the sphere"""
        return self.polygon.centroid

    @property
    def length(self) -> float:
        """diagonal length of the rectangular footprint over the sphere in degrees"""
        return self.polygon.length

    @property
    def circumference(self) -> float:
        """circumference of the rectangular footprint in degrees"""
        return self.polygon.boundary.length

    @property
    def area(self) -> float:
        """area of this footprint on the sphere in degrees squared"""
        return self.polygon.area

    @property
    def possibly_intersecting_projregions(self) -> int:
        """number of possibly intersecting projection regions"""
        if self.area > sc.SkyCell.area:
            return (
                # the number of times the smallest projection region could fit in the image footprint
                np.ceil(self.area / sc.ProjectionRegion.MIN_AREA)
                # plus multiplier for partial intersections on the perimeter
                * 8
            )
        else:
            # 4 foundational intersections
            return 4

    @property
    def possibly_intersecting_skycells(self) -> int:
        """number of possibly intersecting skycells"""
        if self.area > sc.SkyCell.area:
            return (
                # number of times a skycell could fit in the image footprint
                np.ceil(self.area / sc.SkyCell.area)
                # plus multiplier for partial intersections on the perimeter
                * 8
            )
        else:
            # 4 foundational intersections
            return 4

    @property
    def possible_intersecting_projregion_distance(self) -> float:
        """maximum possible distance to the center of an intersecting projection region"""
        return (self.length + sc.ProjectionRegion.MAX_LENGTH) / 2.0

    @property
    def possible_intersecting_skycell_distance(self) -> float:
        """maximum possible distance to the center of an intersecting sky cell"""
        return (self.length + sc.SkyCell.length) / 2.0

    def __str__(self) -> str:
        return f"footprint {self.polygon.boundary.points.to_lonlats()}"

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({self.polygon.boundary.points.to_lonlats()!r})"
        )


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
        sky map instance; defaults to global SKYMAP (Default value = None)

    Returns
    -------
    A sequence of the indices of all sky cells that overlap the supplied image
    (in the referenced sky map). The indices may be used to obtain all
    necessary information about the sky cells, by calling `romancal.skycell.skymap.SkyCell(index)`.
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
        skymap.projection_regions_kdtree.query(
            footprint.center.xyz,
            # NOTE: `k` can be `footprint.possibly_intersecting_projregions` if we assume non-degenerate tesselation
            k=100,
            distance_upper_bound=footprint.possible_intersecting_projregion_distance
            * 1.1,
        )[1]
    )
    nearby_projregion_indices = nearby_projregion_indices[
        np.where(nearby_projregion_indices != len(skymap.model.projection_regions))
    ]

    for projregion_index in nearby_projregion_indices:
        projregion = sc.ProjectionRegion(projregion_index)
        if footprint.polygon.intersects_poly(projregion.polygon):
            # query the LOCAL k-d tree of skycells for possible intersection candidates in (normalized) 3D space
            projregion_nearby_skycell_indices = np.array(
                projregion.skycells_kdtree.query(
                    footprint.center.xyz,
                    # NOTE: `k` can be `footprint.possibly_intersecting_skycells` if we assume non-degenerate tesselation
                    k=100,
                    distance_upper_bound=footprint.possible_intersecting_skycell_distance
                    * 1.1,
                )[1]
            )
            projregion_nearby_skycell_indices = np.array(
                projregion_nearby_skycell_indices[
                    np.where(
                        projregion_nearby_skycell_indices != len(projregion.skycells)
                    )
                ]
            )

            # convert to absolute indices over the full skycell table
            projregion_nearby_skycell_indices += projregion.data["skycell_start"]

            # find polygons that intersect the image footprint
            for skycell_index in projregion_nearby_skycell_indices:
                # print(nearby_skycell_index)
                skycell = sc.SkyCell(skycell_index)
                if footprint.polygon.intersects_poly(skycell.polygon):
                    # print(f"candidate {nearby_skycell_index} intersects")
                    intersecting_skycell_indices.append(skycell_index)

    return np.array(intersecting_skycell_indices)
