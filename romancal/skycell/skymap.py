import logging
import os
import re
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


class SkyCell:
    """
    Square subregion of a projection region, 4.6 arcminutes per side.
    """

    _index: int | None
    _data: np.void
    _skymap: "SkyMap"

    # average area of a sky cell in square degrees on the sphere
    area = 1.7760288493318122e-06

    # average diagonal length of a skycell in degrees on the sphere
    length = 0.0018846729895155984

    def __init__(self, index: int | None, skymap: "SkyMap" = None):
        """
        Parameters
        ----------
        index : int
            Index in the global sky map.
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)
        """

        if skymap is None:
            skymap = SKYMAP

        self._index = index
        self._skymap = skymap
        self._data = None

    @classmethod
    def from_name(cls, name: str, skymap: "SkyMap" = None) -> "SkyCell":
        """
        Retrieve a sky cell from the sky map by its name (see handbook [1] for explanation).

        Parameters
        ----------
        name : str
            Name of a sky cell, for instance `315p86x50y75`.
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)

        References
        ----------
        .. [1] `Skymap Tessellation <https://roman-docs.stsci.edu/data-handbook-home/wfi-data-format/skymap-tessellation>`_
        """
        if not re.match(r"\d{3}\w\d{2}x\d{2}y\d{2}", name):
            raise ValueError(f"invalid skycell name {name}")

        if skymap is None:
            skymap = SKYMAP

        indices = np.where(skymap.model.skycells["name"] == name)[0]
        if len(indices) > 0:
            return SkyCell(indices[0], skymap=skymap)
        else:
            raise KeyError(
                f"skycell with name '{name}' does not exist in the currently-loaded sky map"
            )

    @classmethod
    def from_data(cls, data: np.void, skymap: "SkyMap" = None) -> "SkyCell":
        """
        build an index-less sky cell instance from a data array

        Parameters
        ----------
        data : numpy.void
            array with sky cell parameters (see schema)
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)
        """
        instance = cls(index=None, skymap=skymap)
        instance._data = data
        return instance

    @classmethod
    def from_asn(cls, asn: dict | str, skymap: "SkyMap" = None) -> "SkyCell":
        """
        retrieve a sky cell from WCS info or a target specified in an association

        Attempts to find a sky cell name from the following in order:
            - `skycell_wcs_info.name`
            - `target`

        Parameters
        ----------
        asn : dict | str
            association dictionary or a path to an association file to load
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)
        """

        if isinstance(asn, str | os.PathLike):
            asn = ModelLibrary._load_asn(asn)

        skycell_name = None
        if "skycell_wcs_info" in asn and isinstance(asn["skycell_wcs_info"], dict):
            skycell_name = asn["skycell_wcs_info"]["name"]
        elif "target" in asn and asn["target"].lower() != "none":
            skycell_name = asn["target"]
        else:
            raise ValueError(
                "cannot extract skycell information from an association without `skycell_wcs_info` or `target`"
            )

        return SkyCell.from_name(skycell_name, skymap=skymap)

    @property
    def index(self) -> int | None:
        """index of this sky cell in the sky map"""
        return self._index

    @property
    def data(self) -> np.void:
        """
        Sky cell data.

        ("name", "ra_center", "dec_center", "orientat", "x_tangent", "y_tangent", "ra_corn1", "dec_corn1", "ra_corn2", "dec_corn2", "ra_corn3", "dec_corn3", "ra_corn4", "dec_corn4")
        """
        if self._data is None and self.index is not None:
            self._data = self._skymap.model.skycells[self.index]
        return self._data

    @property
    def name(self) -> str:
        """
        name of this sky cell, for instance `315p86x50y75`

        NOTE
        ----
        the name of a sky cell center comprises the rounded center coordinates of its containing projection region in right ascension and declination,
        and the XY location of the sky cell within its projection region in units of ordinal sky cells from that center
        """
        return self.data["name"].item()

    @property
    def radec_center(self) -> tuple[float, float]:
        """center point in right ascension and declination"""
        return self.data["ra_center"].item(), self.data["dec_center"].item()

    @property
    def orientation(self) -> float:
        return self.data["orientat"].item()

    @property
    def xy_tangent(self) -> tuple[float, float]:
        """center point in pixel coordinates"""
        return self.data["x_tangent"].item(), self.data["y_tangent"].item()

    @property
    def radec_corners(
        self,
    ) -> NDArray[float]:
        """corners in right ascension and declination in the order given by the sky map"""
        return np.array(
            (
                self.data[["ra_corn1", "dec_corn1"]].item(),
                self.data[["ra_corn2", "dec_corn2"]].item(),
                self.data[["ra_corn3", "dec_corn3"]].item(),
                self.data[["ra_corn4", "dec_corn4"]].item(),
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
        return sgv.normalize_vector(np.array(sgv.lonlat_to_vector(*self.radec_center)))

    @cached_property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this sky cell"""
        return sgp.SingleSphericalPolygon(
            points=self.vectorpoint_corners,
            inside=self.vectorpoint_center,
        )

    @cached_property
    def projection_region(self) -> "ProjectionRegion":
        """Projection region containing this sky cell."""
        if self.index is None:
            raise ValueError("no index provided")
        return ProjectionRegion.from_skycell_index(self.index)

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        return self._skymap.pixel_scale

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels across"""
        return self._skymap.pixel_shape

    @property
    def wcs_info(self) -> dict[str, float | str]:
        """WCS properties as defined in the Level 3 association schema"""

        return {
            "name": self.name,
            "pixel_scale": self.pixel_scale,
            "ra_projection_center": self.projection_region.radec_tangent[0],
            "dec_projection_center": self.projection_region.radec_tangent[1],
            "x0_projection": self.xy_tangent[0],
            "y0_projection": self.xy_tangent[1],
            "ra_center": self.radec_center[0],
            "dec_center": self.radec_center[1],
            "nx": self.pixel_shape[0],
            "ny": self.pixel_shape[1],
            "orientat": self.orientation,
            "orientat_projection_center": self.projection_region.orientation
            if self.projection_region.data["dec_max"] != 90.0
            # rotate north polar cap by 180 degrees
            # TODO: find out why this is necessary...
            else self.projection_region.orientation + 180,
        }

    @cached_property
    def wcs(self) -> WCS:
        """WCS representing this sky cell"""
        wcsobj = wcsinfo_to_wcs(
            self.wcs_info,
            bounding_box=(
                (-0.5, self.pixel_shape[0] - 0.5),
                (-0.5, self.pixel_shape[1] - 0.5),
            ),
        )
        wcsobj.array_shape = self.pixel_shape
        return wcsobj

    def __eq__(self, other) -> bool:
        if not isinstance(other, SkyCell):
            return False

        return self.data == other.data


class ProjectionRegion:
    """Projection region in the sky map."""

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
            index of the projection region in the sky map array
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)
        """

        if skymap is None:
            skymap = SKYMAP

        self._index = index
        self._skymap = skymap
        if index is not None:
            self._data = self._skymap.model.projection_regions[index]

    @classmethod
    def from_data(cls, data: np.void, skymap: "SkyMap" = None) -> "ProjectionRegion":
        """
        build a projection region instance from a data array

        Parameters
        ----------
        data : numpy.void
            array with projection region parameters (see schema)
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)
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
            index of the sky cell
        skymap: SkyMap
            sky map instance (defaults to global SKYMAP)

        Returns
        -------
        ProjectionRegion
            projection region corresponding to the given sky cell
        """

        if skymap is None:
            skymap = SKYMAP

        projregion_indices = np.where(
            (index >= skymap.model.projection_regions["skycell_start"])
            & (index < skymap.model.projection_regions["skycell_end"])
        )[0]
        if len(projregion_indices) == 1:
            return cls(projregion_indices[0], skymap=skymap)
        else:
            msg = (
                f"sky cell index {index} not found in any projection regions"
                if len(projregion_indices) == 0
                else f"sky cell index {index} found in multiple projection regions; possibly malformed sky map at {skymap.path}"
            )
            raise KeyError(msg)

    @property
    def index(self) -> int | None:
        """index in the sky map"""
        return self._index

    @property
    def data(self) -> np.void:
        """
        Projection region data.

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
        """indices of sky cells in the sky map within this region"""
        return np.arange(self.data["skycell_start"], self.data["skycell_end"])

    @property
    def skycells(self) -> np.void:
        """subset array of sky cells from the sky map within this region"""
        return self._skymap.model.skycells[self.skycell_indices]

    @cached_property
    def skycells_kdtree(self) -> KDTree:
        """
        LOCAL k-d tree of skycells in this projection region, using normalized center vectorpoints in 3D space

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


class SkyMap:
    """
    Abstract representation of the sky map, comprising of 4058 overlapping rectangular "projection regions" defining gnomonic projection on to uniform pixel grids.
    The pixel scale for all projection regions is identical.

    Each projection region is subdivided into ~2000 square subregions ("sky cells", ~8 million in total), each 4.6' across. These sky cells also overlap each other by a standard number of pixels.

    References
    ----------
    .. [1] `Skymap Tessellation <https://roman-docs.stsci.edu/data-handbook-home/wfi-data-format/skymap-tessellation>`_
    """

    _path: None | Path
    _data: AsdfFile

    def __init__(self, path: None | Path | str = None):
        """
        Parameters
        ----------
        path : None | Path | str, optional
            load sky map from the specified ASDF file (defaults to latest `skycells` ref on CRDS)
        """
        if path is not None and not isinstance(path, Path):
            path = Path(path)
        self._path = path
        self._data = None

    @property
    def path(self) -> None | Path:
        """location of sky map reference file on filesystem"""
        return self._path

    @path.setter
    def path(self, path: None | Path):
        self._path = path
        # reset data if retrieved
        self._data = None

    @property
    def model(self) -> AsdfFile:
        """data model of sky map"""
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
        """
        k-d tree of all skycells in the skymap, using normalized center vectorpoints in 3D space

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
        """number of pixels per sky cell"""
        return self.model.meta.nxy_skycell, self.model.meta.nxy_skycell

    def __getitem__(self, index: int) -> SkyCell:
        """`SkyCell` at the given index in the sky cells array"""
        return SkyCell(index)


SKYMAP = SkyMap(path=os.environ.get("SKYMAP_PATH", None))
