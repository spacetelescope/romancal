"""
Plotting utilities for plotting skycells against supplied image

matplotlib dependency is optional.
"""

import numpy as np
import spherical_geometry.vector as sgv
from numpy.typing import NDArray

import romancal.skycell.skymap as sc

try:
    from matplotlib import pyplot as plt
except ImportError:
    print("matplotlib is required for this plotting utility")

RAD_TO_ARCSEC = 180.0 / np.pi * 3600.0


def find_closest_tangent_point(
    image_corners: list[tuple[float, float]] | tuple[list[float], list[float]],
) -> tuple[tuple[float, float, float], list[int]]:
    """
    Out of all projection regions, find the closest tangent point to the center coordinate of the image.
    """

    # Convert image corners to unit vectors.
    # To deal with the corner case, it is necessary to use spherical_geometry to average the the corners properly.
    image_corner_vectorpoints = sgv.normalize_vector(
        sc.image_coords_to_vec(image_corners)
    )

    # find the closest projection region center to any of the image corners
    corner_nearest = dict(
        sc.SKYMAP.projregions_kdtree.query(corner_vectorpoint)
        for corner_vectorpoint in image_corner_vectorpoints
    )
    closest_projregion_index = corner_nearest[min(corner_nearest)]

    return sc.SKYMAP.projregions[closest_projregion_index][
        ["ra_tangent", "dec_tangent"]
    ], closest_projregion_index


def veccoords_to_tangent_plane(
    vertices: list[tuple[float, float, float]],
    tangent_vectorpoint: tuple[float, float, float],
) -> NDArray[float]:
    """
    Convert the spherical geometry vectors to tangent plane coordinates
    in arcseconds. This algorithm is not precise, but should be good
    enough for now (and besides, the goal here is visualizaion, not
    ultra-precision). This also breaks down numerically very near the
    poles.
    """

    # First compute the tangent plane axis vectors.
    x_axis = sgv.normalize_vector(np.cross([0, 0, 1], tangent_vectorpoint))
    y_axis = sgv.normalize_vector(
        np.array([0, 0, 1])
        - np.array(tangent_vectorpoint)
        * np.dot(np.array([0, 0, 1]), np.array(tangent_vectorpoint))
    )
    avertices = np.vstack(vertices).T
    x_coords = np.dot(x_axis, avertices) * RAD_TO_ARCSEC
    y_coords = np.dot(y_axis, avertices) * RAD_TO_ARCSEC
    return np.stack([x_coords, y_coords], axis=1)


plt.ion()


def plot_field(corners: NDArray[float], id: str = "", fill=None, color=None):
    plt.fill(corners[:, 0], corners[:, 1], color=fill, edgecolor=color)


def plot_skycell(
    skycell: sc.SkyCell,
    tangent_vectorpoint: tuple[float, float, float],
    color=None,
    label: bool = True,
):
    corners_tangentplane = veccoords_to_tangent_plane(
        np.stack(sgv.lonlat_to_vector(*skycell.radec_corners.T), axis=1),
        tangent_vectorpoint,
    )
    corners_tangentplane = np.concatenate(
        [corners_tangentplane, corners_tangentplane[None, -1]]
    )

    plt.plot(corners_tangentplane[:, 0], corners_tangentplane[:, 1], color=color)

    if label:
        center = np.mean(corners_tangentplane[:-1], axis=0)
        plt.annotate(skycell.name, center, va="center", ha="center", size=10)


def plot_image_footprint_and_skycells(
    image_corners: list[tuple[float, float]],
    touched_skycell_indices: list[int],
    nearby_skycell_indices: list[int],
):
    """
    This plots a list of skycell footprints against the image footprint.

    Both the touched skycells as well as nearby skycells are plotted.


    Parameters
    ----------
    image_corners : Either a squence of 4 (ra, dec) pairs, or
        equivalent 2-d numpy array

    touched_skycell_ids: A list of the indices in the skycell table of
        skycells touched by the image footprint.

    candidate_skycell_ids: A list of the indices selected to see if they
        close enough to test if touched by the image footprint.
    """

    plt.clf()
    plt.gca().invert_xaxis()
    plt.plot(0, 0, "*", markersize=10)

    closest_tangent_point, _ = find_closest_tangent_point(image_corners)
    closest_tangent_vectorpoint = sgv.normalize_vector(
        sgv.lonlat_to_vector(*closest_tangent_point)
    )

    image_corners_tangentplane = veccoords_to_tangent_plane(
        sc.image_coords_to_vec(image_corners),
        closest_tangent_vectorpoint,
    )
    plot_field(image_corners_tangentplane, fill="lightgrey", color="black")

    for skycell_index in nearby_skycell_indices:
        plot_skycell(
            sc.SkyCell(skycell_index),
            closest_tangent_vectorpoint,
            color="lightgray",
        )

    for skycell_index in touched_skycell_indices:
        plot_skycell(
            sc.SkyCell(skycell_index),
            closest_tangent_vectorpoint,
            color="blue",
        )

    plt.xlabel("Offset from nearest tangent point in arcsec")
    plt.ylabel("Offset from nearest tangent point in arcsec")

    plt.title(f"{closest_tangent_point} of tangent point in degrees")
