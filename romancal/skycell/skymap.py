import logging
import os
from datetime import datetime
from functools import cached_property
from pathlib import Path

import crds
import numpy as np
import roman_datamodels
import spherical_geometry.great_circle_arc as sga
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from asdf import AsdfFile
from gwcs import WCS
from numpy.typing import NDArray
from scipy.spatial import KDTree

from romancal.datamodels.library import ModelLibrary
from romancal.lib.wcsinfo_to_wcs import wcsinfo_to_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["SKYMAP", "ProjectionRegion", "SkyCells", "SkyMap"]


class SkyCells:
    """Set of square subregions of projection regions, 4.6 arcminutes per side."""

    _indices: NDArray[int]
    _data: np.void
    _skymap: "SkyMap"

    # average area of a skycell in square degrees on the sphere
    area = 1.7760288493318122e-06

    # average diagonal length of a skycell in degrees on the sphere
    length = 0.107984

    def __init__(self, indices: NDArray[int], skymap: "SkyMap" = None):
        """
        Parameters
        ----------
        indices : list[int]
            Indices of skycells in the global skymap.
        skymap: SkyMap
            skymap instance (defaults to global SKYMAP)
        """

        if skymap is None:
            skymap = SKYMAP

        self._indices = np.array(indices)
        self._skymap = skymap
        self._data = None

    @classmethod
    def from_names(cls, names: list[str], skymap: "SkyMap" = None) -> "SkyCells":
        """Retrieve skycells from the sky map [skymap]_ by name.

        Note
        ----
        Order of input names is NOT preserved.

        Parameters
        ----------
        name : list[str]
            List of names of skycells, for instance `315p86x50y75`.
        skymap : SkyMap
            sky map instance; defaults to global SKYMAP (Default value = None)
        """

        if skymap is None:
            skymap = SKYMAP

        if isinstance(names, str):
            names = [names]

        indices = np.isin(skymap.model.skycells["name"], names).nonzero()[0]
        found_names = skymap.model.skycells["name"][indices]

        if len(indices) == len(names):
            return SkyCells(
                indices,
                skymap=skymap,
            )
        else:
            raise KeyError(
                f"no skycells found with the following name(s) in the currently-loaded skymap: {[name for name in names if name not in found_names]}"
            )

    @classmethod
    def from_asns(
        cls, asns: list[os.PathLike | dict], skymap: "SkyMap" = None
    ) -> "SkyCells":
        """retrieve skycells from the given association file(s)

        Attempts to find a skycell name from the following in order:
            - `skycell_wcs_info.name`
            - `target`

        Parameters
        ----------
        asns : list[os.PathLike | dict]
            list of path(s) to association file(s) to load, or dictionaries
        skymap: SkyMap
            skymap instance; defaults to global SKYMAP (Default value = None)
        """

        skycell_names = []
        for asn in asns:
            if isinstance(asn, os.PathLike):
                asn = ModelLibrary._load_asn(asn)

            if "skycell_wcs_info" in asn and isinstance(asn["skycell_wcs_info"], dict):
                skycell_names.append(asn["skycell_wcs_info"]["name"])
            elif "target" in asn and asn["target"].lower() != "none":
                skycell_names.append(asn["target"])
            else:
                raise ValueError(
                    "cannot extract skycell information from an association without `skycell_wcs_info` or `target`"
                )

        return SkyCells.from_names(skycell_names, skymap=skymap)

    @property
    def indices(self) -> NDArray[int]:
        """indices of these skycells in the loaded sky map"""
        return self._indices

    @property
    def data(self) -> np.void:
        """data of these skycells

        ("name", "ra_center", "dec_center", "orientat", "x_tangent", "y_tangent", "ra_corn1", "dec_corn1", "ra_corn2", "dec_corn2", "ra_corn3", "dec_corn3", "ra_corn4", "dec_corn4")
        """
        if self._data is None:
            self._data = self._data = (
                self._skymap.model.skycells[self.indices]
                if len(self.indices) > 0
                else np.array([], dtype=self._skymap.model.skycells.dtype)
            )

        return self._data

    @property
    def names(self) -> list[str]:
        """names of these skycells, for instance `315p86x50y75`

        NOTE
        ----
        the name of a skycell comprises the rounded center coordinates of its containing projection region in right ascension and declination,
        and the XY location of the skycell within its projection region in units of ordinal skycells from that center
        """
        return self.data["name"].tolist()

    @property
    def radec_centers(self) -> NDArray[float]:
        """center points in right ascension and declination (Nx2 array of floats)"""
        return np.array((self.data["ra_center"], self.data["dec_center"])).T

    @property
    def orientations(self) -> NDArray[float]:
        return self.data["orientat"]

    @property
    def xy_tangents(self) -> NDArray[float]:
        """tangent points of projection regions in pixel coordinates (Nx2 array of floats)"""
        return np.array((self.data["x_tangent"], self.data["y_tangent"]))

    @property
    def radec_corners(
        self,
    ) -> NDArray[float]:
        """corners in right ascension and declination in the order given by the loaded sky map (Nx4x2 array of floats)"""
        return np.stack(
            (
                np.stack(
                    [self.data[f"ra_corn{index}"] for index in range(1, 5)], axis=1
                ),
                np.stack(
                    [self.data[f"dec_corn{index}"] for index in range(1, 5)], axis=1
                ),
            ),
            axis=2,
        )

    @cached_property
    def vectorpoint_corners(self) -> NDArray[float]:
        """corners in 3D Cartesian space on the unit sphere (Nx4x3 array of floats)"""
        radec_corners = np.reshape(self.radec_corners, (len(self) * 4, 2))
        return np.reshape(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(radec_corners[:, 0], radec_corners[:, 1]),
                    axis=1,
                )
            ),
            (len(self), 4, 3),
        )

    @cached_property
    def vectorpoint_centers(self) -> NDArray[float]:
        """centers in 3D Cartesian space on the unit sphere (Nx3 array of floats)"""
        return sgv.normalize_vector(
            np.stack(
                sgv.lonlat_to_vector(
                    self.radec_centers[:, 0], self.radec_centers[:, 1]
                ),
                axis=1,
            )
        )

    @cached_property
    def polygons(self) -> sgp.SphericalPolygon:
        """spherical polygons representing these skycells"""
        return sgp.SphericalPolygon(
            [
                sgp.SingleSphericalPolygon(
                    points=vectorpoint_corners,
                    inside=vectorpoint_center,
                )
                for vectorpoint_corners, vectorpoint_center in zip(
                    self.vectorpoint_corners, self.vectorpoint_centers, strict=True
                )
            ]
        )

    @cached_property
    def projection_regions(self) -> list[int]:
        """index of projection region containing each sky cell"""
        skycell_projection_regions = np.empty_like(self.indices)
        for (
            projregion_index,
            projregion_skycell_start_index,
            projregion_skycell_end_index,
        ) in self._skymap.model.projection_regions[
            ["index", "skycell_start", "skycell_end"]
        ]:
            skycell_projection_regions[
                (self.indices >= projregion_skycell_start_index)
                & (self.indices < projregion_skycell_end_index)
            ] = projregion_index
        return skycell_projection_regions

    @property
    def wcs_infos(self) -> list[dict[str, float | str]]:
        """WCS properties as defined in the Level 3 association schema for each skycell"""

        return [
            {
                "name": self._skymap.model.skycells[skycell_index]["name"].item(0),
                "pixel_scale": self.pixel_scale,
                "ra_projection_center": self._skymap.model.projection_regions[
                    "ra_tangent"
                ][projregion_index],
                "dec_projection_center": self._skymap.model.projection_regions[
                    "dec_tangent"
                ][projregion_index],
                "x0_projection": self._skymap.model.skycells["x_tangent"][
                    skycell_index
                ],
                "y0_projection": self._skymap.model.skycells["y_tangent"][
                    skycell_index
                ],
                "ra_center": self._skymap.model.skycells["ra_center"][skycell_index],
                "dec_center": self._skymap.model.skycells["dec_center"][skycell_index],
                "nx": self.pixel_shape[0],
                "ny": self.pixel_shape[1],
                "orientat": self._skymap.model.skycells[skycell_index]["orientat"].item(
                    0
                ),
                "orientat_projection_center": self._skymap.model.projection_regions[
                    "orientat"
                ][projregion_index],
            }
            for skycell_index, projregion_index in zip(
                self.indices, self.projection_regions, strict=True
            )
        ]

    @cached_property
    def wcs(self) -> list[WCS]:
        """WCS objects representing these skycells"""

        wcsobjs = []
        for wcs_info in self.wcs_infos:
            wcsobj = wcsinfo_to_wcs(
                wcs_info,
                bounding_box=(
                    (-0.5, self.pixel_shape[0] - 0.5),
                    (-0.5, self.pixel_shape[1] - 0.5),
                ),
            )
            wcsobj.array_shape = wcsobj.pixel_shape
            wcsobjs.append(wcsobj)
        return wcsobjs

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        return self._skymap.pixel_scale

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels across"""
        return self._skymap.pixel_shape

    def containing(self, radec: NDArray[float]) -> dict[int, list[int]]:
        """
        point(s) contained by each of these skycells

        Parameters
        ----------
        radec: NDArray[float]
            right ascension and declination of coordinate(s)

        Returns
        -------
        mapping of skycell indices to indices of given points contained by that skycell
        """

        radec = np.array(radec)
        if radec.ndim == 1:
            radec = np.expand_dims(radec, axis=0)

        projregion_indices, skycell_projregion_index_indices = np.unique(
            self.projection_regions,
            return_inverse=True,
        )

        projregions: dict[int, list[int]] = {
            projregion_index.item(0): [] for projregion_index in projregion_indices
        }
        for skycell_projregion_index_index, skycell_index in zip(
            skycell_projregion_index_indices, self.indices, strict=True
        ):
            projregions[projregion_indices[skycell_projregion_index_index]].append(
                skycell_index
            )

        skycells: dict[int, list[int]] = {}
        for (
            projregion_index,
            projregion_skycell_indices,
        ) in projregions.items():
            projregion = ProjectionRegion(projregion_index)
            projregion_x, projregion_y = projregion.wcs.world_to_pixel_values(
                radec[:, 0], radec[:, 1]
            )

            for skycell_projregion_index, (
                skycell_x_tangent,
                skycell_y_tangent,
            ) in enumerate(
                self._skymap.model.skycells[projregion_skycell_indices][
                    ["x_tangent", "y_tangent"]
                ]
            ):
                skycell_point_indices = (
                    (
                        projregion_x
                        >= projregion.xy_tangent[0]
                        - skycell_x_tangent
                        - self.pixel_shape[0]
                    )
                    & (
                        projregion_x
                        < projregion.xy_tangent[0]
                        - skycell_x_tangent
                        + self.pixel_shape[0]
                    )
                    & (
                        projregion_y
                        >= projregion.xy_tangent[1]
                        - skycell_y_tangent
                        - self.pixel_shape[1]
                    )
                    & (
                        projregion_y
                        < projregion.xy_tangent[1]
                        - skycell_y_tangent
                        + self.pixel_shape[1]
                    )
                ).nonzero()[0]
                if len(skycell_point_indices) > 0:
                    skycell_index = (
                        skycell_projregion_index + projregion.data["skycell_start"]
                    )
                    if skycell_index not in skycells:
                        skycells[skycell_index] = []
                    skycells[skycell_index].extend(skycell_point_indices)
        return skycells

    def cores_containing(self, radec: NDArray[np.float64]) -> dict[int, list[int]]:
        """
        point(s) exclusively core-contained by each of these skycells

        Parameters
        ----------
        radec: NDArray[float]
            right ascension and declination of coordinate(s)

        Returns
        -------
        mapping of skycell indices to indices of given points exclusively core-contained by that skycell
        """

        radec = np.array(radec)
        if radec.ndim == 1:
            radec = np.expand_dims(radec, axis=0)

        projregion_indices, skycell_projregion_index_indices = np.unique(
            self.projection_regions,
            return_inverse=True,
        )

        projregions: dict[int, list[int]] = {
            projregion_index.item(0): [] for projregion_index in projregion_indices
        }
        for skycell_projregion_index_index, skycell_index in zip(
            skycell_projregion_index_indices, self.indices, strict=True
        ):
            projregions[projregion_indices[skycell_projregion_index_index]].append(
                skycell_index
            )

        skycells: dict[int, list[int]] = {}
        for (
            projregion_index,
            projregion_skycell_indices,
        ) in projregions.items():
            projregion = ProjectionRegion(projregion_index)
            projregion_points_within = projregion.contains_radec(radec)
            # only continue if any points lie within the projection region
            if np.any(projregion_points_within):
                projregion_radec = radec[projregion_points_within]
                projregion_x, projregion_y = projregion.wcs.invert(
                    projregion_radec[:, 0],
                    projregion_radec[:, 1],
                )

                projregion_skycells = SkyCells(projregion_skycell_indices)
                for projregion_skycell_index, (
                    skycell_name,
                    skycell_x_tangent,
                    skycell_y_tangent,
                ) in zip(
                    projregion_skycells.indices,
                    projregion_skycells.data[["name", "x_tangent", "y_tangent"]],
                    strict=True,
                ):
                    skycell_x = projregion_x - (
                        projregion.data["x_tangent"] - skycell_x_tangent
                    )
                    skycell_y = projregion_y - (
                        projregion.data["y_tangent"] - skycell_y_tangent
                    )

                    half_margin = self._skymap.model.meta["skycell_border_pixels"] / 2
                    core_contains = (
                        (half_margin - 0.5 < skycell_x)
                        & (skycell_x < self.pixel_shape[0] - half_margin - 0.5)
                        & (half_margin - 0.5 < skycell_y)
                        & (skycell_y < self.pixel_shape[1] - half_margin - 0.5)
                    )

                    # handle polar singularities
                    if np.any(np.abs(radec[projregion_points_within, 1]) == 90):
                        # TODO if the polar projection regions change, this will need to be updated
                        if str(skycell_name).endswith("x50y50"):
                            core_contains[
                                np.abs(radec[projregion_points_within, 1]) == 90
                            ] = True

                    if np.any(core_contains):
                        projregion_skycell_index = projregion_skycell_index.item(0)
                        if projregion_skycell_index not in skycells:
                            skycells[projregion_skycell_index] = []
                        skycells[projregion_skycell_index].extend(
                            value.item(0)
                            for value in projregion_points_within.nonzero()[0][
                                core_contains
                            ]
                        )
        return skycells

    @cached_property
    def kdtree(self) -> KDTree:
        """k-d tree of skycells, using normalized center vectorpoints in 3D space"""
        return KDTree(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(
                        self.data["ra_center"],
                        self.data["dec_center"],
                    ),
                    axis=1,
                )
            )
        )

    def __len__(self) -> int:
        return len(self._indices)

    def __eq__(self, other) -> bool:
        if isinstance(other, SkyCells):
            return self.indices == other.indices
        else:
            return False

    def __str__(self) -> str:
        return f"{self.names} [{self._skymap.path}]"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.indices}, {self._skymap})"


class ProjectionRegion:
    """Projection region in the loaded sky map."""

    _index: int | None
    _data: np.void
    _skymap: "SkyMap"

    # area of the smallest projection region in square degrees on the sphere
    #   min(sc.ProjectionRegion(index).polygon.area() for index in range(len(sc.SKYMAP.model.projection_regions)))
    MIN_AREA = 0.002791388883915502

    # diagonal length of the longest projection region in degrees on the sphere
    #   max(sc.ProjectionRegion(index).length for index in range(len(sc.SKYMAP.model.projection_regions)))
    MAX_LENGTH = 0.08174916691321586

    def __init__(self, index: int | None, skymap: "SkyMap" = None):
        """
        Parameters
        ----------
        index : int
            index of the projection region in the loaded sky map array
        skymap: SkyMap
            skymap instance (defaults to global SKYMAP)
        """

        if skymap is None:
            skymap = SKYMAP

        self._index = index
        self._skymap = skymap
        if index is not None:
            self._data = self._skymap.model.projection_regions[index]

    @classmethod
    def from_data(cls, data: np.void, skymap: "SkyMap" = None) -> "ProjectionRegion":
        """build a projection region instance from a data array

        Parameters
        ----------
        data : numpy.void
            array with projection region parameters (see schema)
        skymap: SkyMap
            skymap instance; defaults to global SKYMAP (Default value = None)
        """
        instance = cls(index=data["index"], skymap=skymap)
        instance._data = data
        return instance

    @classmethod
    def from_skycell_index(
        cls, index: int, skymap: "SkyMap" = None
    ) -> "ProjectionRegion":
        """
        Parameters
        ----------
        index : int
            index of the skycell
        skymap : SkyMap
            skymap instance; defaults to global SKYMAP (Default value = None)

        Returns
        -------
        ProjectionRegion
            projection region corresponding to the given skycell
        """

        if skymap is None:
            skymap = SKYMAP

        projregion_indices = (
            (index >= skymap.model.projection_regions["skycell_start"])
            & (index < skymap.model.projection_regions["skycell_end"])
        ).nonzero()[0]
        if len(projregion_indices) == 1:
            return cls(projregion_indices[0], skymap=skymap)
        else:
            msg = (
                f"skycell index {index} not found in any projection regions"
                if len(projregion_indices) == 0
                else f"skycell index {index} found in multiple projection regions; possibly malformed sky map at {skymap.path}"
            )
            raise KeyError(msg)

    @property
    def index(self) -> int | None:
        """index in the loaded sky map"""
        return self._index

    @property
    def data(self) -> np.void:
        """Projection region data.

        ("index", "ra_tangent", "dec_tangent", "ra_min", "ra_max", "dec_min", "dec_max", "orientat", "x_tangent", "y_tangent", "nx", "ny", "skycell_start", "skycell_end")
        """
        return self._data

    @property
    def radec_tangent(self) -> tuple[float, float]:
        """projection origin (tangent point with the celestial sphere) in right ascension and declination"""
        return self.data[["ra_tangent", "dec_tangent"]].item()

    @property
    def radec_bounds(self) -> tuple[float, float, float, float]:
        """bounds in right ascension and declination in order [xmin, ymin, xmax, ymax]"""
        return self.data[["ra_min", "dec_min", "ra_max", "dec_max"]].item()

    @property
    def orientation(self) -> float:
        return self.data["orientat"].item()

    @property
    def xy_tangent(self) -> tuple[float, float]:
        """projection origin (tangent point with the celestial sphere) in pixel coordinates"""
        return self.data[["x_tangent", "y_tangent"]].item()

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels across"""
        return self.data[["nx", "ny"]].item()

    @cached_property
    def skycell_indices(self) -> NDArray[int]:
        """indices of sky cells in the loaded sky map within this projection region"""
        return np.arange(self.data["skycell_start"], self.data["skycell_end"])

    @cached_property
    def skycells(self) -> SkyCells:
        """collection of all skycells in this projection region"""
        return SkyCells(self.skycell_indices)

    @property
    def radec_corners(
        self,
    ) -> NDArray:
        """corners in right ascension and declination in clockwise order"""
        return np.array(
            (
                self.data[["ra_min", "dec_min"]].item(),
                self.data[["ra_max", "dec_min"]].item(),
                self.data[["ra_max", "dec_max"]].item(),
                self.data[["ra_min", "dec_max"]].item(),
            )
        )

    @cached_property
    def vectorpoint_corners(self) -> NDArray[float]:
        """corners in 3D Cartesian space on the unit sphere"""
        return sgv.normalize_vector(
            np.stack(sgv.lonlat_to_vector(*np.array(self.radec_corners).T), axis=1)
        )

    @cached_property
    def vectorpoint_center(self) -> tuple[float, float, float]:
        """center in 3D Cartesian space on the unit sphere"""
        return np.mean(self.vectorpoint_corners, axis=0)

    @cached_property
    def length(self) -> float:
        """diagonal length of the region"""
        # assume radial against sky background
        return max(
            sga.length(
                self.vectorpoint_corners[index], self.vectorpoint_corners[index + 2]
            )
            for index in range(len(self.vectorpoint_corners) - 3)
        )

    @property
    def is_polar(self) -> bool:
        """whether this projection region is a polar cap"""
        return self.data["dec_max"] == 90.0 or self.data["dec_min"] == -90.0

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this region"""
        if self.is_polar:
            # the projection regions at the poles are circular caps on the sphere;
            # a polygon built from the corners in that case would be degenerate
            return sgp.SingleSphericalPolygon.from_cone(
                *self.radec_tangent,
                radius=(
                    90.0 - self.radec_bounds[1]
                    if self.data["dec_max"] == 90.0
                    else -self.radec_bounds[3]
                ),
                steps=16,
            )
        else:
            return sgp.SingleSphericalPolygon(
                points=self.vectorpoint_corners,
                inside=self.vectorpoint_center,
            )

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        return self._skymap.pixel_scale

    @property
    def wcs_info(self) -> dict[str, float | str]:
        """WCS properties as defined in the Level 3 association schema"""

        return {
            "name": f"proj{self.index}",
            "pixel_scale": self.pixel_scale,
            "ra_projection_center": self.data["ra_tangent"],
            "dec_projection_center": self.data["dec_tangent"],
            "x0_projection": self.data["x_tangent"],
            "y0_projection": self.data["y_tangent"],
            "nx": self.data["nx"],
            "ny": self.data["ny"],
            "orientat": self.orientation,
        }

    @cached_property
    def wcs(self) -> WCS:
        """WCS representing this skycell"""
        wcsobj = wcsinfo_to_wcs(
            self.wcs_info,
            bounding_box=(
                (-0.5, self.pixel_shape[0] - 0.5),
                (-0.5, self.pixel_shape[1] - 0.5),
            ),
        )
        wcsobj.array_shape = self.pixel_shape
        return wcsobj

    def contains_radec(self, radec: NDArray[float]) -> NDArray[bool]:
        """whether the given point(s) are contained within the bounds of this projection region"""
        radec = np.array(radec)
        if radec.ndim == 1:
            radec = np.expand_dims(radec, axis=0)

        # only one projection region contains each point
        return (
            ra_in_range(
                radec[:, 0],
                self.data["ra_min"],
                self.data["ra_max"],
            )
            & (radec[:, 1] >= self.data["dec_min"])
            & (radec[:, 1] < self.data["dec_max"])
        )

    def __eq__(self, other) -> bool:
        if not isinstance(other, ProjectionRegion):
            return False

        return self.data == other.data

    def __str__(self) -> str:
        return f"projregion {self.index} [{self._skymap.path}]"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.index}, {self._skymap})"


class SkyMap:
    """Abstract representation of the skymap [skymap]_, comprising of 4058 overlapping rectangular
    "projection regions" defining gnomonic projection on to uniform pixel grids.
    The pixel scale for all projection regions is identical.

    Each projection region is subdivided into ~2000 square subregions ("skycells", ~8 million in total),
    each 4.6' across. These skycells also overlap each other by a standard number of pixels.

    References
    ----------
    .. [skymap] `Skymap Tessellation <https://roman-docs.stsci.edu/data-handbook-home/wfi-data-format/skymap-tessellation>`_
    """

    _path: None | Path
    _data: AsdfFile

    def __init__(self, path: None | Path | str = None):
        """
        Parameters
        ----------
        path : None | Path | str, optional
            load skymap from the specified ASDF file (defaults to latest `skycells` ref on CRDS)
        """
        if path is not None and not isinstance(path, Path):
            path = Path(path)
        self._path = path
        self._data = None

    @property
    def path(self) -> None | Path:
        """location of skymap reference file on filesystem"""
        return self._path

    @path.setter
    def path(self, path: None | Path):
        self._path = path
        # reset data if retrieved
        self._data = None

    @property
    def model(self) -> AsdfFile:
        """data model of skymap"""
        if self._data is None:
            if self._path is None:
                rmap = crds.getreferences(
                    {
                        "roman.meta.instrument.name": "WFI",
                        "roman.meta.exposure.start_time": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
                    },
                    reftypes=["skycells"],
                    observatory="roman",
                )
                self._path = Path(rmap["skycells"])
            self._data = roman_datamodels.open(self._path, memmap=True)
        return self._data

    @cached_property
    def skycells(self) -> SkyCells:
        """collection of all skycells in this skymap"""
        return SkyCells(np.arange(len(self.model.skycells)))

    @cached_property
    def projection_regions_kdtree(self) -> KDTree:
        """k-d tree of all projection regions in this skymap, using normalized center vectorpoints in 3D space"""
        return KDTree(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(
                        self.model.projection_regions["ra_tangent"],
                        self.model.projection_regions["dec_tangent"],
                    ),
                    axis=1,
                )
            )
        )

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        # The scale is given in arcseconds per pixel. Convert to degrees.
        return self.model.meta["pixel_scale"] / 3600.0

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels per skycell"""
        return self.model.meta.nxy_skycell, self.model.meta.nxy_skycell

    def projection_regions_containing(
        self, radec: NDArray[float]
    ) -> dict[int, list[int]]:
        """
        point(s) contained by each projection region in this sky map

        Parameters
        ----------
        radec: NDArray[float]
            right ascension and declination of coordinate(s)

        Returns
        -------
        mapping of projection region indices to indices of given points contained by that projection region
        """

        radec = np.array(radec)
        if radec.ndim == 1:
            radec = np.expand_dims(radec, axis=0)

        projregions = {}
        for (
            projregion_index,
            ra_min,
            ra_max,
            dec_min,
            dec_max,
        ) in self.model.projection_regions[
            ["index", "ra_min", "ra_max", "dec_min", "dec_max"]
        ]:
            # each point SHOULD only be contained by a single projection region
            contained_point_indices = (
                ra_in_range(
                    radec[:, 0],
                    ra_min,
                    ra_max,
                )
                & (radec[:, 1] >= dec_min)
                & (radec[:, 1] < dec_max)
            ).nonzero()[0]
            if len(contained_point_indices) > 0:
                projregions[projregion_index] = contained_point_indices

        return projregions

    def __getitem__(self, indices: int) -> SkyCells:
        """`SkyCells` at the given indices in the sky cells array"""
        return SkyCells(indices)

    def __str__(self) -> str:
        return f"skymap {self.path}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.path})"


def ra_in_range(ra: float, low: float, high: float):
    """whether the given longitude lies within the given min and max range, handling wrapping"""
    ra = ra % 360
    low = low % 360
    high = high % 360
    if high == low:
        high = 360.0
    if low <= high:
        return (ra >= low) & (ra <= high)
    else:
        return (ra >= low) | (ra <= high)


SKYMAP = SkyMap(path=os.environ.get("SKYMAP_PATH", None))
