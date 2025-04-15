"""
Plotting utilities for plotting skycells against supplied image

matplotlib dependency is optional.
"""

import numpy as np
import spherical_geometry.vector as sgv
from numpy.typing import NDArray
from spherical_geometry.vector import normalize_vector

from .skycells import SKYCELLS, ProjectionRegion, image_coords_to_vec

try:
    from matplotlib import pyplot as plt
except ImportError:
    print("matplotlib is required for this plotting utility")

RAD_TO_ARCSEC = 180.0 / np.pi * 3600.0


def get_cartesian_corners(
    projection_region: ProjectionRegion,
) -> tuple[NDArray[float], NDArray[float]]:
    """
    Construct vertex coordinates for a projection region definition suitable
    for plotting the defined region (ra coordinates, dec coordinates). It returns
    coordinates in the cartesian system so that it avoids the distortions in plots
    due to areas approaching the poles.
    """

    corners = np.array(projection_region.radec_corners)
    return sgv.lonlat_to_vector(corners[:, 0], corners[:, 1])


def find_closest_tangent_point(
    projection_regions: list[ProjectionRegion],
    image_corners: list[tuple[float, float]] | tuple[list[float], list[float]],
) -> tuple[tuple[float, float, float], list[int]]:
    """
    Out of all listed projection regions, find the closest tangent point to the center coordinate of the image.
    """

    # To deal with the corner case, it is necessary to use spherical_geometry to average the the corners properly.

    # Convert image corners to unit vectors.
    image_corner_vectorpoints = image_coords_to_vec(image_corners)
    image_center_vectorpoint = normalize_vector(image_corner_vectorpoints.mean(axis=0))

    unique_tangent_points = []
    for projection_region in projection_regions:
        tangent_point = projection_region.radec_center
        if tangent_point not in unique_tangent_points:
            unique_tangent_points.append(tangent_point)

    # Compute 3D cartesian distance for each tangent point from the image center
    distances = [
        (
            (image_center_vectorpoint - np.array(sgv.lonlat_to_vector(*tangent_point)))
            ** 2
        ).sum()
        for tangent_point in unique_tangent_points
    ]
    sorted_distance_indices = sorted(
        zip(distances, range(len(distances)), strict=False)
    )

    # find closest tangent point by distance
    closest_tangent_point = np.array(
        sgv.lonlat_to_vector(*unique_tangent_points[sorted_distance_indices[0][1]])
    )

    # Now associate sorted distance indices with those of all projection regions
    projregion_tangentpoint_indices = [
        index
        for projection_region in projection_regions
        for index, (_, tangent_point_index) in enumerate(sorted_distance_indices)
        if unique_tangent_points[tangent_point_index] == projection_region.radec_center
    ]
    return closest_tangent_point, projregion_tangentpoint_indices


def veccoords_to_tangent_plane(
    vertices: list[tuple[float, float, float]],
    tangent_point_vec: tuple[float, float, float],
) -> tuple[list[float], list[float]]:
    """
    Convert the spherical geometry vectors to tangent plane coordinates
    in arcseconds. This algorithm is not precise, but should be good
    enough for now (and besides, the goal here is visualizaion, not
    ultra-precision). This also breaks down numerically very near the
    poles.
    """
    # First compute the tangent plane axis vectors.
    x_axis = normalize_vector(np.cross([0, 0, 1], tangent_point_vec))
    y_axis = normalize_vector(
        np.array([0, 0, 1])
        - np.array(tangent_point_vec)
        * np.dot(np.array([0, 0, 1]), np.array(tangent_point_vec))
    )
    avertices = np.vstack(vertices)
    x_coords = np.dot(x_axis, avertices) * RAD_TO_ARCSEC
    y_coords = np.dot(y_axis, avertices) * RAD_TO_ARCSEC
    return x_coords, y_coords


plt.ion()


def plot_field(corners, id="", fill=None, color=None):
    plt.fill(corners[0], corners[1], color=fill, edgecolor=color)


def plot_skycell(corners, id="", color=None):
    plt.plot(corners[0], corners[1], color=color)
    if id:
        idstr = str(SKYCELLS.skycells[id]["index"])
        center = (corners[0][:-1].mean(), corners[1][:-1].mean())
        plt.annotate(idstr, center, va="center", ha="center", size=10)


def plot(image_corners, touched_skycell_indices, candidate_skycell_indices):
    """
    This plots a list of skycell footprints against the image footprint.

    Both the touched skycells as well as candidate skycells are plotted.


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
    touched_skycells = [SKYCELLS.skycells[index] for index in touched_skycell_indices]
    candidate_skycells = [
        SKYCELLS.skycells[index] for index in candidate_skycell_indices
    ]
    tangent_point, touched_skycell_index = find_closest_tangent_point(
        touched_skycells, image_corners
    )
    ra, dec = sgv.vector_to_lonlat(*tangent_point)
    dummy, skycell_tp_index = find_closest_tangent_point(
        candidate_skycells, image_corners
    )
    vec_image_corners = image_coords_to_vec(image_corners)
    tp_image_corners = veccoords_to_tangent_plane(vec_image_corners, tangent_point)
    plot_field(tp_image_corners, fill="lightgrey", color="black")
    for candidate_skycell, candidate_skycell_index in zip(
        candidate_skycells, candidate_skycell_indices, strict=False
    ):
        plot_skycell(
            veccoords_to_tangent_plane(
                get_cartesian_corners(candidate_skycell), tangent_point
            ),
            id=candidate_skycell_index,
            color="lightgray",
        )
    for touched_skycell, touched_skycell_index in zip(
        touched_skycells, touched_skycell_indices, strict=False
    ):
        plot_skycell(
            veccoords_to_tangent_plane(
                get_cartesian_corners(touched_skycell), tangent_point
            ),
            id=touched_skycell_index,
            color="blue",
        )
    plt.xlabel("Offset from nearest tangent point in arcsec")
    plt.ylabel("Offset from nearest tangent point in arcsec")
    plt.title(f"RA: {ra} Dec: {dec} of tangent point in degrees")
