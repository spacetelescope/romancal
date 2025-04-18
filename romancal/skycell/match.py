"""
This module determines which projection regions (associated with "sky tiles") overlap with the given image.

Currently this assumes that the sky projected borders of all images are straight.
"""

import logging
from functools import cached_property

import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from gwcs import WCS
from numpy.typing import NDArray

from .skymap import SKYMAP, SkyCell, SkyTile, image_coords_to_vec

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Radius in degrees within which projection region centers must lie within to be tested for intersection
CANDIDATE_RADIUS = 0.5


class ImageFootprint:
    __radec_corners: NDArray[float]

    def __init__(
        self,
        *radec_corners: tuple[float, float],
    ):
        if len(radec_corners) != 4:
            raise ValueError(f"need 4 corners, not {len(radec_corners)}")
        self.__radec_corners = np.array(radec_corners)

    @classmethod
    def from_gwcs(
        cls, iwcs: WCS, image_shape: tuple[int, int] | None = None
    ) -> "ImageFootprint":
        # Now must find size of corresponding image, with three possible sources of that information.
        if image_shape is None:
            # Both bounding_box and pixel_shape are in x, y order contrary to numpy convention.
            if hasattr(iwcs, "bounding_box") and iwcs.bounding_box is not None:
                # Presumes that the bounding_box matches image array boundaries
                bbintervals = iwcs.bounding_box.intervals
                # This compensates for the half pixel adjustment in the general code.
                image_shape = (bbintervals[1].upper + 0.5, bbintervals[0].upper + 0.5)
            elif hasattr(iwcs, "pixel_shape") and iwcs.pixel_shape is not None:
                image_shape = (iwcs.pixel_shape[1], iwcs.pixel_shape[0])
            else:
                raise ValueError(
                    "Use of a wcs object requires at least one of the bounding_box"
                    " or pixel_shape attributes be set to the image shape or that the"
                    "image_shape argument be set"
                )

        # Compute the image corners ra, dec from the wcs
        (cxm, cxp), (cym, cyp) = (
            (-0.5, image_shape[1] - 0.5),
            (-0.5, image_shape[0] - 0.5),
        )

        return cls(iwcs(cxp, cyp), iwcs(cxm, cyp), iwcs(cxm, cym), iwcs(cxp, cym))

    @property
    def radec_corners(
        self,
    ) -> NDArray:
        return self.__radec_corners

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        corner_vectorpoints = image_coords_to_vec(self.radec_corners)

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
        footprint = ImageFootprint(*image_corners)

    # TODO: There are 8 million sky cells in the table, ~2000 per skytile.
    # TODO: We should try to derive the maximum number of intersecting skycells `n` based on the image size,
    # TODO: and then use a spatial index to find the `n` nearest skycells to the image center to test intersection with the image footprint.
    nearby_skycell_indices = []
    intersecting_skycell_indices = []
    for skytile_index in range(len(SKYMAP.skytiles)):
        skytile = SkyTile(skytile_index)
        if footprint.polygon.intersects_poly(skytile.polygon):
            skycells = SKYMAP.skycells[skytile.skycell_indices]

            # Convert all celestial coordinates to cartesion coordinates.
            skycell_center_vectorpoints = np.array(
                sgv.lonlat_to_vector(
                    skycells[:]["ra_center"], skycells[:]["dec_center"]
                )
            ).transpose()
            image_corner_vectorpoints = image_coords_to_vec(footprint.radec_corners)

            # Approximate center of image by averaging corner vectors
            image_center_vectorpoint = sgv.normalize_vector(
                np.mean(image_corner_vectorpoints, axis=0)
            )

            # Compute distances in 3D space between image center and projection region centers
            skycell_center_distances = np.sqrt(
                np.sum(
                    (skycell_center_vectorpoints - image_center_vectorpoint) ** 2,
                    axis=1,
                )
            )

            # find skycell indices below candidate radius
            skytile_nearby_skycell_indices = np.where(
                skycell_center_distances < CANDIDATE_RADIUS
            )
            nearby_skycell_indices.extend(skytile_nearby_skycell_indices[0])

            # find polygons that intersect the image footprint
            for skycell_candidate_index in skytile_nearby_skycell_indices[0]:
                # make index absolute to entire table
                skycell_candidate_index += skytile.data["skycell_start"]
                # print(skycell_candidate_index)
                skycell = SkyCell(skycell_candidate_index)
                if footprint.polygon.intersects_poly(skycell.polygon):
                    # print(f"candidate {skycell_candidate_index} intersects")
                    intersecting_skycell_indices.append(skycell_candidate_index)

    return intersecting_skycell_indices, nearby_skycell_indices
