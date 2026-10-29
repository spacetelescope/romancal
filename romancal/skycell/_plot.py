"""
Plotting utilities for plotting skycells against supplied image

matplotlib dependency is optional.
"""

import numpy as np
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
    "plot_skycells",
    "radec_to_tangent_plane",
]

DEG_TO_ARCSEC = 3600.0


def find_intersecting_projregions(
    footprint: sm._ImageFootprint,
    skymap: sc.SkyMap = None,
) -> list[int]:
    """Out of all projection regions, find ones with skycells that intersect the given image footprint

    Parameters
    ----------
    footprint: sm._ImageFootprint :
        sequence of points (ra, dec) or an `_ImageFootprint` object
    skymap: sc.SkyMap :
        skymap instance; defaults to global SKYMAP (Default value = None)

    Returns
    -------
    indices of projection regions in the skymap that intersect the given footprint
    """

    if skymap is None:
        skymap = sc.SKYMAP

    skycells = sc.SkyCells(
        np.array(
            sm.find_skycell_matches(footprint.radec_corners, skymap=skymap), dtype=int
        ),
        skymap=skymap,
    )
    return np.unique(skycells.projection_regions).tolist()


def radec_to_tangent_plane(
    radec: NDArray[float], radec_tangent: tuple[float, float]
) -> NDArray[float]:
    """Project (..., 2) right ascension and declination onto the plane tangent
    to the sky at `radec_tangent`, as (..., 2) offsets in arcseconds.

    This is the gnomonic projection used by the skymap, so skycell edges and
    other great circles are straight lines.
    """
    radec = np.asarray(radec)
    x, y = sm._tangent_plane(*radec_tangent)(radec[..., 0], radec[..., 1])
    return np.stack((x, y), axis=-1) * DEG_TO_ARCSEC


def _closed(points: NDArray[float]) -> NDArray[float]:
    """repeat the first of the given (..., N, 2) points at the end"""
    return np.concatenate([points, points[..., :1, :]], axis=-2)


def _per_skycell(values, count: int) -> list:
    """one value per skycell, from nothing, a single value, or one per skycell"""
    if values is None or isinstance(values, str):
        return [values] * count
    if len(values) == 1:
        return list(values) * count
    return list(values)


def plot_field(corners: NDArray[float], id: str = "", fill=None, color=None, axis=None):
    if axis is None:
        axis = plt
    axis.fill(corners[:, 0], corners[:, 1], color=fill, edgecolor=color)


def plot_projregion(
    projregion: sc.ProjectionRegion, color=None, label: bool = True, axis=None
):
    """plot the nominal right ascension and declination bounds of a projection region"""
    if axis is None:
        axis = plt

    ra_min, dec_min, ra_max, dec_max = projregion.radec_bounds
    if ra_max <= ra_min:
        ra_max += 360
    ra = np.linspace(ra_min, ra_max, 100)
    dec = np.linspace(dec_min, dec_max, 100)
    if projregion.is_polar:
        # a polar cap is bounded by a single circle of declination
        cap_dec = dec_min if dec_max == 90 else dec_max
        boundary = np.stack((ra, np.full_like(ra, cap_dec)), axis=-1)
    else:
        # lines of constant declination are curved in the tangent plane
        boundary = np.concatenate(
            [
                np.stack((ra, np.full_like(ra, dec_min)), axis=-1),
                np.stack((np.full_like(dec, ra_max), dec), axis=-1),
                np.stack((ra[::-1], np.full_like(ra, dec_max)), axis=-1),
                np.stack((np.full_like(dec, ra_min), dec[::-1]), axis=-1),
            ]
        )
    boundary = radec_to_tangent_plane(_closed(boundary), projregion.radec_tangent)

    axis.plot(boundary[:, 0], boundary[:, 1], color=color)

    if label:
        axis.annotate(
            f"proj{projregion.index}",
            (0, 0),
            va="center",
            ha="center",
            size=10,
            color=color,
        )


def plot_skycells(
    skycells: sc.SkyCells,
    radec_tangent: tuple[float, float],
    colors=None,
    labels: list[str] | None = None,
    annotations: list[str] | None = None,
    axis=None,
):
    """plot the outlines of skycells on the plane tangent at `radec_tangent`

    `colors` and `labels` may be a single value, or one per skycell.
    """
    if axis is None:
        axis = plt

    colors = _per_skycell(colors, len(skycells))
    labels = _per_skycell(labels, len(skycells))
    outlines = radec_to_tangent_plane(_closed(skycells.radec_corners), radec_tangent)

    for index, outline in enumerate(outlines):
        axis.plot(
            outline[:, 0], outline[:, 1], color=colors[index], label=labels[index]
        )

        if annotations:
            axis.annotate(
                annotations[index],
                np.mean(outline[:-1], axis=0),
                va="center",
                ha="center",
                size=10,
                color=colors[index],
            )


def plot_image_footprint_and_skycells(
    footprint: list[tuple[float, float]] | sm._ImageFootprint,
    skycells: sc.SkyCells,
    skymap: sc.SkyMap = None,
) -> list[tuple[Axis, tuple[float, float]]]:
    """This plots a list of skycell footprints against the image footprint.

    Both the touched skycells as well as nearby skycells are plotted, on the
    tangent plane of each projection region with skycells that intersect the
    footprint.

    Parameters
    ----------
    footprint : list | sm._ImageFootprint :
        sequence of points (ra, dec) or an `_ImageFootprint` object
    skycells : sc.SkyCells :
        skycells to highlight, usually those that intersect the footprint
    skymap : sc.SkyMap :
        skymap instance; defaults to global SKYMAP (Default value = None)

    Returns
    -------
    the axis for each projection region, with its tangent point
    """

    if not isinstance(footprint, sm._ImageFootprint):
        footprint = sm._ImageFootprint(footprint)

    if skymap is None:
        skymap = sc.SKYMAP

    # plot each intersecting projection region's intersection onto a tangent plane
    intersecting_projregion_indices = find_intersecting_projregions(
        footprint, skymap=skymap
    )
    axes = []
    for projregion_index in intersecting_projregion_indices:
        projregion = sc.ProjectionRegion(projregion_index, skymap=skymap)
        radec_tangent = projregion.radec_tangent

        figure, axis = plt.subplots(1, 1)
        figure.suptitle(f"projection region {projregion_index}")
        # east to the left, as on the sky
        axis.invert_xaxis()
        axis.set_aspect("equal")
        axis.plot(0, 0, "+", markersize=10)

        plot_projregion(projregion, color="tab:blue", axis=axis)

        plot_skycells(projregion.skycells, radec_tangent, colors="darkgrey", axis=axis)

        projregion_intersecting_skycells = sc.SkyCells(
            skycells.indices[skycells.projection_regions == projregion_index],
            skymap=skymap,
        )
        plot_skycells(
            projregion_intersecting_skycells,
            radec_tangent,
            colors="red",
            annotations=projregion_intersecting_skycells.names,
            axis=axis,
        )

        plot_field(
            radec_to_tangent_plane(footprint.radec_corners, radec_tangent),
            fill="lightgrey",
            color="black",
            axis=axis,
        )

        axis.set_xlabel("Offset from tangent point in arcsec")
        axis.set_ylabel("Offset from tangent point in arcsec")

        axis.set_title(f"tangent point radec {np.array(radec_tangent)}")

        axes.append((axis, radec_tangent))

    return axes
