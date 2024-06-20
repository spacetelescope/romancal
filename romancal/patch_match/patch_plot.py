"""
Plotting utilities for plotting patches against supplied image

matplotlib dependency is optional.
"""

import spherical_geometry.vector as sgv

import romancal.patch_match.patch_match as pm

try:
    from matplotlib import pyplot as plt
except ImportError:
    print("matplotlib is required for this plotting utility")

plt.ion()


def plot_field(corners, id="", fill=None, color=None):
    plt.fill(corners[0], corners[1], color=fill, edgecolor=color)


def plot_patch(corners, id="", color=None):
    plt.plot(corners[0], corners[1], color=color)
    if id:
        idstr = str(pm.PATCH_TABLE[id]["index"])
        center = (corners[0][:-1].mean(), corners[1][:-1].mean())
        plt.annotate(idstr, center, va="center", ha="center", size=10)


def plot(image_corners, patches_touched_ids, patches_candidate_ids):
    """
    This plots a list of patch footprints against the image footprint.

    Both the touched patchs as well as candidate patches are plotted.


    Parameters
    ----------
    image_corners : Either a squence of 4 (ra, dec) pairs, or
        equivalent 2-d numpy array

    patches_touched_ids: A list of the indices into the patch table of
        patches touched by the image footprint.

    patches_candidate_ids: A list of the indices selected to see if they
        close enough to test if touched by the image footprint.
    """
    plt.clf()
    plt.gca().invert_xaxis()
    plt.plot(0, 0, "*", markersize=10)
    patches_touched = [pm.PATCH_TABLE[index] for index in patches_touched_ids]
    patches_candidate = [pm.PATCH_TABLE[index] for index in patches_candidate_ids]
    tangent_point, patch_tp_id_touched = pm.find_closest_tangent_point(
        patches_touched, image_corners
    )
    ra, dec = sgv.vector_to_lonlat(*tangent_point)
    dummy, patch_tp_id = pm.find_closest_tangent_point(patches_candidate, image_corners)
    vec_image_corners = pm.image_coords_to_vec(image_corners)
    tp_image_corners = pm.veccoords_to_tangent_plane(vec_image_corners, tangent_point)
    plot_field(tp_image_corners, fill="lightgrey", color="black")
    for patch, id in zip(patches_candidate, patches_candidate_ids):
        plot_patch(
            pm.veccoords_to_tangent_plane(
                pm.get_cartesian_corners(patch), tangent_point
            ),
            id=id,
            color="lightgray",
        )
    for patch, id in zip(patches_touched, patches_touched_ids):
        plot_patch(
            pm.veccoords_to_tangent_plane(
                pm.get_cartesian_corners(patch), tangent_point
            ),
            id=id,
            color="blue",
        )
    plt.xlabel("Offset from nearest tangent point in arcsec")
    plt.ylabel("Offset from nearest tangent point in arcsec")
    plt.title(f"RA: {ra} Dec: {dec} of tangent point in degrees")
