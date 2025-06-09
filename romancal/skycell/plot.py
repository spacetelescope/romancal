"""
Plotting utilities for plotting skycells against supplied image

matplotlib dependency is optional.
"""

import numpy as np
import spherical_geometry.vector as sgv
from numpy.typing import NDArray

import romancal.skycell.match as sm
import romancal.skycell.skymap as sc

try:
    from matplotlib import pyplot as plt
    from matplotlib.axis import Axis
except ImportError:
    print("matplotlib is required for this plotting utility")

__all__ = [
    "plot_field",
    "plot_image_footprint_and_skycells",
    "plot_projregion",
    "plot_skycell",
    "veccoords_to_tangent_plane",
]

RAD_TO_ARCSEC = 180.0 / np.pi * 3600.0


def find_intersecting_projregions(
    footprint: sm.ImageFootprint,
    skymap: sc.SkyMap = None,
) -> list[int]:
    """Out of all projection regions, find ones that intersect the given image footprint

    Parameters
    ----------
    footprint: sm.ImageFootprint :
        sequence of 4 points (ra, dec) or an `ImageFootprint` object
    skymap: sc.SkyMap :
        sky map instance; defaults to global SKYMAP (Default value = None)

    Returns
    -------
    indices of projection regions in the sky map that intersect the given footprint
    """

    if skymap is None:
        skymap = sc.SKYMAP

    # find the closest projection regions to the image center
    _, nearby_projregion_indices = skymap.projection_regions_kdtree.query(
        footprint.vectorpoint_center,
        k=footprint.possibly_intersecting_projregions,
        distance_upper_bound=footprint.possible_intersecting_projregion_distance * 1.1,
    )

    intersecting_projregion_indices = []
    for projregion_index in nearby_projregion_indices:
        if projregion_index < skymap.projection_regions_kdtree.n:
            projregion = sc.ProjectionRegion(projregion_index, skymap=skymap)
            if footprint.polygon.intersects_poly(projregion.polygon):
                intersecting_projregion_indices.append(projregion_index)

    return intersecting_projregion_indices


def veccoords_to_tangent_plane(
    vertices: list[tuple[float, float, float]],
    tangent_vectorpoint: tuple[float, float, float],
) -> NDArray[float]:
    """Convert the spherical geometry vectors to tangent plane coordinates
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


def plot_field(corners: NDArray[float], id: str = "", fill=None, color=None, axis=None):
    if axis is None:
        axis = plt
    axis.fill(corners[:, 0], corners[:, 1], color=fill, edgecolor=color)


def plot_projregion(
    projregion: sc.ProjectionRegion, color=None, label: bool = True, axis=None
):
    if axis is None:
        axis = plt

    tangent_vectorpoint = sgv.normalize_vector(
        sgv.lonlat_to_vector(*projregion.radec_tangent)
    )
    corners = projregion.vectorpoint_corners
    corners = np.concat([corners, corners[0, :].reshape((1, 3))], axis=0)
    corners_tangentplane = veccoords_to_tangent_plane(
        corners,
        tangent_vectorpoint,
    )

    axis.plot(corners_tangentplane[:, 0], corners_tangentplane[:, 1], color=color)

    if label:
        center = np.mean(corners_tangentplane[:-1], axis=0)
        axis.annotate(
            f"proj{projregion.index}",
            center,
            va="center",
            ha="center",
            size=10,
            color=color,
        )


def plot_skycell(
    skycell: sc.SkyCell,
    tangent_vectorpoint: tuple[float, float, float],
    color=None,
    label: str | None = None,
    annotation: str | None = None,
    axis=None,
):
    if axis is None:
        axis = plt

    corners = skycell.vectorpoint_corners
    corners = np.concat([corners, corners[0, :].reshape((1, 3))], axis=0)
    corners_tangentplane = veccoords_to_tangent_plane(
        corners,
        tangent_vectorpoint,
    )

    axis.plot(
        corners_tangentplane[:, 0], corners_tangentplane[:, 1], color=color, label=label
    )

    if annotation:
        center = np.mean(corners_tangentplane[:-1], axis=0)
        axis.annotate(
            annotation, center, va="center", ha="center", size=10, color=color
        )


def plot_image_footprint_and_skycells(
    footprint: list[tuple[float, float]] | sm.ImageFootprint,
    skymap: sc.SkyMap = None,
) -> list[tuple[Axis, tuple[float, float, float]]]:
    """This plots a list of skycell footprints against the image footprint.

    Both the touched skycells as well as nearby skycells are plotted.

    Parameters
    ----------
    footprint : list | sm.ImageFootprint :
        sequence of 4 points (ra, dec) or an `ImageFootprint` object
    skymap : sc.SkyMap :
        sky map instance; defaults to global SKYMAP (Default value = None)
    """

    if not isinstance(footprint, sm.ImageFootprint):
        footprint = sm.ImageFootprint(footprint)

    if skymap is None:
        skymap = sc.SKYMAP

    # plot each intersecting projection region's intersection onto a tangent plane
    intersecting_projregion_indices = find_intersecting_projregions(
        footprint, skymap=skymap
    )
    axes = []
    for projregion_index in intersecting_projregion_indices:
        figure = plt.figure()
        figure.gca().invert_xaxis()
        figure.suptitle(f"projection region {projregion_index}")
        axis = figure.subplots(1, 1)
        axis.plot(0, 0, "+", markersize=10)

        projregion = sc.ProjectionRegion(projregion_index, skymap=skymap)

        tangent_vectorpoint = sgv.normalize_vector(
            sgv.lonlat_to_vector(*projregion.radec_tangent)
        )
        image_corners_tangentplane = veccoords_to_tangent_plane(
            footprint.vectorpoint_vertices,
            tangent_vectorpoint,
        )
        plot_field(image_corners_tangentplane, fill="lightgrey", color="black")

        plot_projregion(projregion, color="lightgrey")

        for skycell_index in projregion.skycell_indices:
            skycell = sc.SkyCell(skycell_index, skymap=skymap)
            plot_skycell(
                skycell,
                tangent_vectorpoint,
                color="darkgrey",
            )

        _, nearby_skycell_indices = projregion.skycells_kdtree.query(
            footprint.vectorpoint_center,
            k=footprint.possibly_intersecting_skycells,
            distance_upper_bound=footprint.length + sc.SkyCell.length,
        )

        for skycell_index in nearby_skycell_indices:
            skycell = sc.SkyCell(
                projregion.data["skycell_start"] + skycell_index, skymap=skymap
            )
            plot_skycell(
                skycell,
                tangent_vectorpoint,
                color="red",
            )

            if footprint.polygon.intersects_poly(skycell.polygon):
                plot_skycell(
                    skycell, tangent_vectorpoint, color="blue", annotation=skycell.name
                )

        axis.set_xlabel("Offset from nearest tangent point in arcsec")
        axis.set_ylabel("Offset from nearest tangent point in arcsec")

        axis.set_title(f"tangent point radec {np.array(projregion.radec_tangent)}")

        axes.append((axis, tangent_vectorpoint))

    return axes
