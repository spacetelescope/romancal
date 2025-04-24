"""
This module determines which sky cells overlap with the given image.

Currently this assumes that the sky projected borders of all calibrated L2 images are straight.
"""

import logging
from functools import cached_property

import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from gwcs import WCS
from numpy.typing import NDArray

import romancal.skycell.skymap as sm
from romancal.assign_wcs.utils import wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ImageFootprint:
    """abstraction of an image footprint"""

    __radec_bounds: tuple[float, float, float, float]

    def __init__(
        self,
        ra_min: float,
        dec_min: float,
        ra_max: float,
        dec_max: float,
    ):
        """
        Parameters
        ----------
        ra_min : float
            minimum right ascension
        dec_min : float
            minimum declination
        ra_max : float
            maximum right ascension
        dec_max : float
            maximum declination
        """
        self.__radec_bounds = ra_min, dec_min, ra_max, dec_max

    @classmethod
    def from_corners(cls, radec_corners: list[tuple[float, float]]) -> "ImageFootprint":
        """build footprint from corner points in right ascension and declination"""
        radec_corners = np.array(radec_corners)
        if radec_corners.shape[0] != 4:
            raise ValueError(f"need 4 corners, not {radec_corners.shape[0]}")
        return cls(*np.min(radec_corners, axis=0), *np.max(radec_corners, axis=0))

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

        return cls.from_corners(iwcs.footprint(center=False))

    @property
    def radec_bounds(self) -> tuple[float, float, float, float]:
        """bounds in right ascension and declination in order [xmin, ymin, xmax, ymax]"""
        return self.__radec_bounds

    @property
    def radec_corners(self) -> NDArray:
        """corners in right ascension and declination in counterclockwise order"""
        return np.array(
            [
                (self.radec_bounds[0], self.radec_bounds[1]),
                (self.radec_bounds[2], self.radec_bounds[1]),
                (self.radec_bounds[2], self.radec_bounds[3]),
                (self.radec_bounds[0], self.radec_bounds[3]),
            ]
        )

    @cached_property
    def radec_center(self) -> tuple[float, float]:
        """center point in right ascension and declination"""
        return sgv.vector_to_lonlat(
            *sm.image_coords_to_vec(self.radec_corners).mean(axis=0)
        )

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this image footprint"""

        corner_vectorpoints = sm.image_coords_to_vec(self.radec_corners)

        # Approximate center of image by averaging corner vectors
        center_vectorpoint = corner_vectorpoints.mean(axis=0)

        # construct polygon from corner points and center point
        return sgp.SingleSphericalPolygon(
            points=sgv.normalize_vector(corner_vectorpoints),
            inside=sgv.normalize_vector(center_vectorpoint),
        )


def find_skycell_matches(
    image_corners: list[tuple[float, float]] | NDArray[float] | WCS,
    image_shape: tuple[int, int] | None = None,
) -> tuple[list[int], list[int]]:
    """Find projection regions that the image overlaps with

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

    Returns
    -------
    A sequence of the indices of all projection regions that overlap the supplied image
    (in the referenced projection region table). The indices may be used to obtain all
    necessary information about the projection regions.
    """

    if isinstance(image_corners, WCS):
        footprint = ImageFootprint.from_gwcs(image_corners, image_shape)
    else:
        footprint = ImageFootprint.from_corners(image_corners)

    footprint_center_vectorpoint = sgv.normalize_vector(
        sm.image_coords_to_vec(footprint.radec_corners).mean(axis=0)
    )

    # derive the maximum number of possibly intersecting skytiles (analogous to projection regions) and skycells
    # based on the image footprint area and the known skytile and skycell area
    # TODO: improve these calculations, they are probably over-estimating...
    max_num_intersecting_skytiles = 8 + 2 * int(
        np.ceil(footprint.polygon.area() / sm.SKYTILE_AREA)
    )
    max_num_intersecting_skycells = 8 + 2 * int(
        np.ceil(footprint.polygon.area() / sm.SKYCELL_AREA)
    )

    nearby_skycell_indices = []
    intersecting_skycell_indices = []

    # query the global k-d tree of projection regions for possible intersection candidates in (normalized) 3D space
    _, nearby_projregion_indices = sm.SKYMAP.projregions_kdtree.query(
        footprint_center_vectorpoint, k=max_num_intersecting_skytiles
    )

    for projregion_index in nearby_projregion_indices:
        projregion = sm.ProjectionRegion(projregion_index)
        if footprint.polygon.intersects_poly(projregion.polygon):
            # query the LOCAL k-d tree of skycells for possible intersection candidates in (normalized) 3D space
            _, projregion_nearby_skycell_indices = projregion.skycells_kdtree.query(
                footprint_center_vectorpoint, k=max_num_intersecting_skycells
            )

            # convert to absolute indices over the full skycell table
            projregion_nearby_skycell_indices = [
                projregion.data["skycell_start"] + index
                for index in projregion_nearby_skycell_indices
            ]
            nearby_skycell_indices.extend(projregion_nearby_skycell_indices)

            # find polygons that intersect the image footprint
            for nearby_skycell_index in projregion_nearby_skycell_indices:
                # print(nearby_skycell_index)
                skycell = sm.SkyCell(nearby_skycell_index)
                if footprint.polygon.intersects_poly(skycell.polygon):
                    # print(f"candidate {nearby_skycell_index} intersects")
                    intersecting_skycell_indices.append(nearby_skycell_index)

    return intersecting_skycell_indices, nearby_skycell_indices
