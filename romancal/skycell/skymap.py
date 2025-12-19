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
    length = 0.0018846729895155984

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
        """Retrieve skycells from the skymap [skymap]_ by name.

        Note
        ----
        Order of input names is NOT preserved.

        Parameters
        ----------
        name : list[str]
            List of names of skycells, for instance `315p86x50y75`.
        skymap : SkyMap
            skymap instance; defaults to global SKYMAP (Default value = None)

        References
        ----------
        .. [skymap] `Skymap Tessellation <https://roman-docs.stsci.edu/data-handbook-home/wfi-data-format/skymap-tessellation>`_
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
        """indices of these skycells in the skymap"""
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
        """corners in right ascension and declination in the order given by the skymap (Nx4x2 array of floats)"""
        return np.reshape(
            (
                np.concatenate([self.data[f"ra_corn{index}"] for index in range(1, 5)]),
                np.concatenate(
                    [self.data[f"dec_corn{index}"] for index in range(1, 5)]
                ),
            ),
            (2, 4, len(self)),
        ).T

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
        """index of each skycell's projection region"""

        return [
            (
                (
                    skycell_index
                    >= self._skymap.model.projection_regions["skycell_start"]
                )
                & (skycell_index < self._skymap.model.projection_regions["skycell_end"])
            ).nonzero()[0][0]
            for skycell_index in self.indices
        ]

    @property
    def wcs_infos(self) -> list[dict[str, float | str]]:
        """WCS properties as defined in the Level 3 association schema for each skycell"""

        return [
            {
                "name": self._skymap.model.skycells[skycell_index]["name"].item(0),
                "pixel_scale": self.pixel_scale,
                "ra_projection_center": self._skymap.model.projection_regions[
                    skycell_projection_region
                ]["ra_tangent"].item(0),
                "dec_projection_center": self._skymap.model.projection_regions[
                    skycell_projection_region
                ]["dec_tangent"].item(0),
                "x0_projection": self._skymap.model.skycells[skycell_index][
                    "x_tangent"
                ].item(0),
                "y0_projection": self._skymap.model.skycells[skycell_index][
                    "y_tangent"
                ].item(0),
                "ra_center": self._skymap.model.skycells[skycell_index][
                    "ra_center"
                ].item(0),
                "dec_center": self._skymap.model.skycells[skycell_index][
                    "dec_center"
                ].item(0),
                "nx": self.pixel_shape[0],
                "ny": self.pixel_shape[1],
                "orientat": self._skymap.model.skycells[skycell_index]["orientat"].item(
                    0
                ),
                "orientat_projection_center": self._skymap.model.projection_regions[
                    skycell_projection_region
                ]["orientat"].item(0),
            }
            for skycell_index, skycell_projection_region in zip(
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

    def __getitem__(self, indices: int) -> "SkyCells":
        """`SkyCells` at the given indices in the skycells array"""
        return SkyCells(self.indices[np.isin(self.indices, indices)])

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
    """Projection region in the skymap."""

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
            index of the projection region in the skymap array
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
        """index in the skymap"""
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
        """indices of skycells in the sky map within this region"""
        return np.arange(self.data["skycell_start"], self.data["skycell_end"])

    @property
    def skycells(self) -> np.void:
        """subset array of skycells from the sky map within this region"""
        return self._skymap.model.skycells[self.skycell_indices]

    @cached_property
    def skycells_kdtree(self) -> KDTree:
        """LOCAL k-d tree of skycells in this projection region, using normalized center vectorpoints in 3D space

        NOTE
        ----
        add `skycell_start` to the indices returned by this tree to convert to skycell indices in the parent skymap
        """

        return KDTree(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(
                        self.skycells["ra_center"],
                        self.skycells["dec_center"],
                    ),
                    axis=1,
                )
            )
        )

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
    def skycells_kdtree(self) -> KDTree:
        """k-d tree of all skycells in the skymap, using normalized center vectorpoints in 3D space

        NOTE
        ----
        there are 8 million skycells in the skymap; constructing this tree will take a long time. It is recommended that you instead use the `.skycells_kdtree` property of an individual projection region instead.
        """
        return KDTree(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(
                        self.model.skycells["ra_center"],
                        self.model.skycells["dec_center"],
                    ),
                    axis=1,
                )
            )
        )

    @cached_property
    def projection_regions_kdtree(self) -> KDTree:
        """k-d tree of all projection regions in the skymap, using normalized center vectorpoints in 3D space"""
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

    def __getitem__(self, indices: int) -> SkyCells:
        """`SkyCells` at the given indices in the skycells array"""
        return SkyCells(indices)

    def __str__(self) -> str:
        return f"skymap {self.path}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.path})"


SKYMAP = SkyMap(path=os.environ.get("SKYMAP_PATH", None))
