"""
Plotting utilities for plotting skycells against supplied image

matplotlib dependency is optional.
"""

import spherical_geometry.vector as sgv

from .match import (
    SKYCELLS_TABLE,
    find_closest_tangent_point,
    get_cartesian_corners,
    image_coords_to_vec,
    veccoords_to_tangent_plane,
)

try:
    from matplotlib import pyplot as plt
except ImportError:
    print("matplotlib is required for this plotting utility")

plt.ion()


def plot_field(corners, id="", fill=None, color=None):
    plt.fill(corners[0], corners[1], color=fill, edgecolor=color)


def plot_skycell(corners, id="", color=None):
    plt.plot(corners[0], corners[1], color=color)
    if id:
        idstr = str(SKYCELLS_TABLE["roman"]["skycells"][id]["index"])
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
    touched_skycells = [
        SKYCELLS_TABLE["roman"]["skycells"][index] for index in touched_skycell_indices
    ]
    candidate_skycells = [
        SKYCELLS_TABLE["roman"]["skycells"][index]
        for index in candidate_skycell_indices
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
