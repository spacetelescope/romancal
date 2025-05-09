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
from romancal.assign_wcs.utils import wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ImageFootprint:
    """abstraction of an image footprint"""

    _radec_corners: NDArray[float]

    def __init__(self, radec_corners: list[tuple[float, float]]):
        """
        Parameters
        ----------
        radec_corners: list[tuple[float, float]]
            four corners of the image in right ascension and declination
        """
        radec_corners = np.array(radec_corners)
        if radec_corners.shape[0] != 4:
            raise ValueError(f"need 4 corners, not {radec_corners.shape[0]}")
        self._radec_corners = radec_corners

    @classmethod
    def from_gwcs(
        cls, iwcs: WCS, image_shape: tuple[int, int] | None = None
    ) -> "ImageFootprint":
        """create an image footprint from a GWCS object (and image shape, if no bounding box is present)"""

        if iwcs.bounding_box is None:
            if image_shape is None:
                # Now must find size of corresponding image, with three possible sources of that information.
                # Both bounding_box and pixel_shape are in x, y order contrary to numpy convention.
                if hasattr(iwcs, "bounding_box") and iwcs.bounding_box is not None:
                    bbintervals = iwcs.bounding_box.intervals
                    # This compensates for the half pixel adjustment in the general code.
                    image_shape = (
                        bbintervals[1].upper + 0.5,
                        bbintervals[0].upper + 0.5,
                    )
                elif hasattr(iwcs, "pixel_shape") and iwcs.pixel_shape is not None:
                    image_shape = (iwcs.pixel_shape[1], iwcs.pixel_shape[0])
                else:
                    raise ValueError(
                        "Cannot infer image footprint from GWCS object because "
                        "`image_shape` was not specified and "
                        "the GWCS object does not have `.bounding_box` nor `.pixel_shape`."
                    )

            iwcs.bounding_box = wcs_bbox_from_shape(image_shape)

        return cls(iwcs.footprint(center=False))

    @property
    def radec_corners(self) -> NDArray:
        """corners in right ascension and declination in counterclockwise order"""
        return self._radec_corners

    @cached_property
    def radec_center(self) -> tuple[float, float]:
        """center point in right ascension and declination"""
        return sgv.vector_to_lonlat(*self.vectorpoint_center)

    @cached_property
    def vectorpoint_corners(self) -> NDArray[float]:
        """corners in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(
            np.stack(sgv.lonlat_to_vector(*np.array(self.radec_corners).T), axis=1)
        )

    @cached_property
    def vectorpoint_center(self) -> tuple[float, float, float]:
        """center in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(np.mean(self.vectorpoint_corners, axis=0))

    @cached_property
    def length(self) -> float:
        """diagonal length of the rectangular footprint"""
        # assume radial against sky background
        return sga.length(self.vectorpoint_corners[0], self.vectorpoint_corners[2])

    @cached_property
    def circumference(self) -> float:
        """circumference of the rectangular footprint"""
        return sum(
            sga.length(
                self.vectorpoint_corners[index], self.vectorpoint_corners[index + 1]
            )
            for index in range(-1, len(self.vectorpoint_corners) - 1)
        )

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this image footprint"""
        return sgp.SingleSphericalPolygon(
            points=self.vectorpoint_corners,
            inside=self.vectorpoint_center,
        )

    @cached_property
    def area(self) -> float:
        """area of this footprint on the sphere in degrees squared"""
        return self.polygon.area()

    @cached_property
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

    @cached_property
    def possibly_intersecting_skycells(self) -> int:
        """number of possibly intersecting skycells"""
        if self.polygon.area() > sc.SkyCell.area:
            return (
                # number of times a skycell could fit in the image footprint
                np.ceil(self.polygon.area() / sc.SkyCell.area)
                # plus multiplier for partial intersections on the perimeter
                * 8
            )
        else:
            # 4 foundational intersections
            return 4


def find_skycell_matches(
    image_corners: list[tuple[float, float]] | NDArray[float] | WCS,
    image_shape: tuple[int, int] | None = None,
    skymap: sc.SkyMap = None,
) -> tuple[list[int], list[int]]:
    """
    Find sky cells overlapping the provided image footprint

    Parameters
    ----------
    image_corners : Either a squence of 4 (ra, dec) pairs, or
        equivalent 2-d numpy array, or
        a GWCS instance. The instance must have either the bounding_box or
        pixel_shape attribute defined, or the following image_shape argument
        must be supplied
    image_shape : image shape to be used if a GWCS instance is supplied
        and does not contain a value for either the bounding_box or
        pixel_shape attributes. Default value is None.
    skymap: SkyMap
        sky map instance (defaults to global SKYMAP)

    Returns
    -------
    A sequence of the indices of all sky cells that overlap the supplied image
    (in the referenced sky map). The indices may be used to obtain all
    necessary information about the sky cells, by calling `romancal.skycell.skymap.SkyCell(index)`.
    """

    if isinstance(image_corners, WCS):
        footprint = ImageFootprint.from_gwcs(image_corners, image_shape)
    else:
        footprint = ImageFootprint(image_corners)

    if skymap is None:
        skymap = sc.SKYMAP

    intersecting_skycell_indices = []

    # query the global k-d tree of projection regions for possible intersection candidates in (normalized) 3D space
    nearby_projregion_indices = np.array(
        skymap.projection_regions_kdtree.query(
            footprint.vectorpoint_center,
            k=footprint.possibly_intersecting_projregions,
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
                    footprint.vectorpoint_center,
                    k=footprint.possibly_intersecting_skycells,
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
