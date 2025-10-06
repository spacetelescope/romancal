import numpy as np
import spherical_geometry.vector as sgv
from matplotlib import pyplot as plt

from romancal.skycell.plot import (
    plot_projregion,
    plot_skycell,
    veccoords_to_tangent_plane,
)
from romancal.skycell.skymap import ProjectionRegion, SkyCell

if __name__ == "__main__":
    skycell_names = [
        "000p86x27y29",
        "000p86x28y30",
        "000p86x29y31",
        "000p86x30y32",
        "000p86x31y32",
        "000p86x32y32",
        "000p86x30y33",
        "000p86x31y33",
        "000p86x32y33",
        "000p86x30y34",
        "000p86x31y34",
        "000p86x32y34",
        "000p86x53y68",
        "000p86x54y68",
        "000p86x55y68",
        "000p86x53y69",
        "000p86x54y69",
        "000p86x60y60",
        "000p86x64y64",
        "135p90x26y48",
        "135p90x25y49",
        "135p90x26y49",
        "135p90x25y50",
        "135p90x26y50",
    ]

    projregion_axes = {}
    for skycell_name in skycell_names:
        try:
            skycell = SkyCell.from_name(skycell_name)
        except KeyError:
            continue

        projregion = skycell.projection_region

        if projregion.index not in projregion_axes:
            projregion_axes[projregion.index] = plt.figure().subplots(1, 1)

            plt.tight_layout()
        projregion_axis = projregion_axes[projregion.index]

        tangent_vectorpoint = sgv.normalize_vector(
            sgv.lonlat_to_vector(*projregion.radec_tangent)
        )

        corners = skycell.vectorpoint_corners
        corners = np.concat([corners, corners[0, :].reshape((1, 3))], axis=0)
        corners_tangentplane = veccoords_to_tangent_plane(
            corners,
            tangent_vectorpoint,
        )

        projregion_axis.imshow(
            skycell.core[::-1, :],
            extent=[
                np.min(corners_tangentplane[:, 0]),
                np.max(corners_tangentplane[:, 0]),
                np.min(corners_tangentplane[:, 1]),
                np.max(corners_tangentplane[:, 1]),
            ],
            alpha=np.where(skycell.core[::-1, :], 1, 0).astype(np.float64),
        )

    for projregion_index, projregion_axis in projregion_axes.items():
        projregion = ProjectionRegion(projregion_index)
        tangent_vectorpoint = sgv.normalize_vector(
            sgv.lonlat_to_vector(*projregion.radec_tangent)
        )
        for skycell_index in projregion.skycell_indices:
            skycell = SkyCell(skycell_index)
            plot_skycell(
                skycell,
                tangent_vectorpoint,
                color="lightgrey",
                axis=projregion_axis,
            )
        plot_projregion(projregion, color="black", axis=projregion_axis)

    plt.show()
