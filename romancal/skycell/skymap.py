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
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames
from numpy.typing import NDArray
from roman_datamodels import stnode
from scipy.spatial import KDTree
from stcal.alignment import util as wcs_util

from romancal.datamodels.library import ModelLibrary

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class SkyCell:
    """
    Square subregion of a projection region, 4.6 arcminutes per side.
    """

    _index: int | None
    _data: np.void

    # average area of a sky cell in square degrees on the sphere
    area = 1.7760288493318122e-06

    # average diagonal length of a skycell in degrees on the sphere
    length = 0.0018846729895155984

    def __init__(self, index: int | None):
        """
        Parameters
        ----------
        index : int
            Index in the global sky map.
        """
        self._index = index
        self._data = None

    @classmethod
    def from_name(cls, name: str) -> "SkyCell":
        """
        Retrieve a sky cell from the sky map by its name (see handbook [1] for explanation).

        Parameters
        ----------
        name : str
            Name of a sky cell, for instance `315p86x50y75`.

        References
        ----------
        .. [1] `Skymap Tessellation <https://roman-docs.stsci.edu/data-handbook-home/wfi-data-format/skymap-tessellation>`_
        """
        if not re.match(r"\d{3}\w\d{2}x\d{2}y\d{2}", name):
            raise ValueError(f"invalid skycell name {name}")

        indices = np.where(SKYMAP.skycells["name"] == name)[0]
        if len(indices) > 0:
            return SkyCell(indices[0])
        else:
            raise KeyError(
                f"skycell with name '{name}' does not exist in the currently-loaded sky map"
            )

    @classmethod
    def from_data(cls, data: np.void) -> "SkyCell":
        """
        build an index-less sky cell instance from a data array

        Parameters
        ----------
        data : numpy.void
            array with sky cell parameters (see schema)
        """
        instance = cls(index=None)
        instance._data = data
        return instance

    @classmethod
    def from_center_and_coordinates(
        cls,
        skytile_center: tuple[float, float],
        skycell_coordinates: tuple[int, int],
    ) -> "SkyCell":
        """
        Retrieve a sky cell from the sky map using
        - the center of its containing sky tile
        - the ordinal XY coordinates of the sky cell from the sky tile origin (number of sky cells away from the center)

        (see handbook [1] for further explanation)

        Parameters
        ----------
        projregion_center : tuple[float, float]
            center coordinates of its containing sky tile in right ascension and declination
        skycell_coordinates : tuple[int, int]
            XY location of the sky cell within its sky tile, in units of ordinal sky cells from the center

        References
        ----------
        .. [1] `Skymap Tessellation <https://roman-docs.stsci.edu/data-handbook-home/wfi-data-format/skymap-tessellation>`_
        """
        return cls.from_name(
            f"r{round(skytile_center[0]):03}d{'p' if skytile_center[1] >= 0 else 'm'}{round(skytile_center[1]):02}x{'p' if skycell_coordinates[1] >= 0 else 'm'}{skycell_coordinates[1]:02}y{'p' if skycell_coordinates[1] >= 0 else 'm'}{skycell_coordinates[1]:02}"
        )

    @classmethod
    def from_asn(cls, asn: dict | Path) -> "SkyCell":
        """
        retrieve a sky cell from WCS info or a target specified in an association

        Parameters
        ----------
        asn : asdf.AsdfTree | Path
            association dictionary or a path to an association file to load
        """
        if isinstance(asn, Path):
            asn = ModelLibrary._load_asn(asn).asn
        if "skycell_wcs_info" in asn and asn["skycell_wcs_info"] != "none":
            skycell_name = asn["skycell_wcs_info"]["name"]
        elif "target" in asn:
            # check to see if the product name contains a skycell name & if true get the skycell record
            skycell_name = asn["target"]
        else:
            raise ValueError(
                "cannot extract skycell information from modellibrary association with neither WCS nor target info"
            )

        return SkyCell.from_name(skycell_name)

    @property
    def index(self) -> int | None:
        """index of this sky cell in the sky map"""
        return self._index

    @property
    def data(self) -> np.void:
        """properties from the sky map"""
        if self._data is None and self.index is not None:
            self._data = SKYMAP.skycells[self.index]
        return self._data

    @property
    def name(self) -> str:
        """
        name of this sky cell, for instance `315p86x50y75`

        NOTE
        ----
        the name of a sky cell center comrpises the rounded coordinates of its containing sky tile in right ascension and declination,
        and the XY location of the sky cell within its sky tile in units of ordinal sky cells
        """
        return self.data[0]

    @property
    def radec_center(self) -> tuple[float, float]:
        """center point in right ascension and declination"""
        return self.data["ra_center"], self.data["dec_center"]

    @property
    def orientation(self) -> float:
        return self.data["orientat"]

    @property
    def xy_tangent(self) -> tuple[float, float]:
        """center point in pixel coordinates"""
        return self.data["x_tangent"], self.data["y_tangent"]

    @property
    def radec_corners(
        self,
    ) -> NDArray[float]:
        """corners in right ascension and declination in the order given by the sky map"""
        return np.array(
            (
                (self.data["ra_corn1"], self.data["dec_corn1"]),
                (self.data["ra_corn2"], self.data["dec_corn2"]),
                (self.data["ra_corn3"], self.data["dec_corn3"]),
                (self.data["ra_corn4"], self.data["dec_corn4"]),
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
        """projection region containing this sky cell"""
        if self.index is None:
            raise ValueError("no index provided")
        return ProjectionRegion.from_skycell_index(self.index)

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        return SKYMAP.pixel_scale

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels across"""
        return SKYMAP.pixel_shape

    @property
    def wcsinfo(self) -> dict:
        """WCS properties"""

        orientation = self.projection_region.orientation
        if self.projection_region.is_polar:
            # rotate polar cap by 180 degrees
            # TODO: find out why this is necessary...
            orientation += 180

        return {
            "ra_ref": self.projection_region.data["ra_tangent"],
            "dec_ref": self.projection_region.data["dec_tangent"],
            "x_ref": self.data["x_tangent"],
            "y_ref": self.data["y_tangent"],
            "rotation_matrix": None,
            "orientat": orientation,
            "pixel_scale": self.pixel_scale,
            "pixel_shape": self.pixel_shape,
            "ra_center": self.data["ra_center"],
            "dec_center": self.data["dec_center"],
            "ra_corn1": self.data["ra_corn1"],
            "dec_corn1": self.data["dec_corn1"],
            "ra_corn2": self.data["ra_corn2"],
            "dec_corn2": self.data["dec_corn2"],
            "ra_corn3": self.data["ra_corn3"],
            "dec_corn3": self.data["dec_corn3"],
            "ra_corn4": self.data["ra_corn4"],
            "dec_corn4": self.data["dec_corn4"],
        }

    @cached_property
    def wcs(self) -> WCS:
        """WCS representing this sky cell"""
        return wcsinfo_to_wcs(self.wcsinfo)

    def __eq__(self, other) -> bool:
        if not isinstance(other, SkyCell):
            return NotImplemented

        return self.data == other.data


class ProjectionRegion:
    """projection region in the sky map, corresponding to a single sky tile on the sky map"""

    _index: int | None
    _data: np.void

    # area of the smallest projection region in square degrees on the sphere
    MIN_AREA = 0.002791388883915502

    # diagonal length of the longest projection region in degrees on the sphere
    MAX_LENGTH = 0.08174916691321586

    def __init__(self, index: int | None):
        """
        Parameters
        ----------
        index : int
            index of the projection region in the sky map array
        """
        self._index = index
        if index is not None:
            self._data = SKYMAP.projection_regions[index]

    @classmethod
    def from_data(cls, data: np.void) -> "ProjectionRegion":
        """
        build an index-less projection region instance from a data array

        Parameters
        ----------
        data : numpy.void
            array with projection region parameters (see schema)
        """
        instance = cls(index=None)
        instance._data = data
        return instance

    @classmethod
    def from_skycell_index(cls, index: int) -> "ProjectionRegion":
        """
        Parameters
        ----------
        index : int
            index of the sky cell

        Returns
        -------
        ProjectionRegion
            projection region corresponding to the given sky cell
        """
        for projregion in SKYMAP.projection_regions:
            if (
                int(projregion["skycell_start"])
                <= index
                < int(projregion["skycell_end"])
            ):
                return cls(projregion["index"])
        else:
            raise KeyError(
                f"sky cell index {index} not found in any projection regions"
            )

    @property
    def index(self) -> int | None:
        """index in the sky map"""
        return self._index

    @property
    def data(self) -> np.void:
        """properties from the sky map"""
        return self._data

    @property
    def radec_tangent(self) -> tuple[float, float]:
        """projection origin (tangent point with the celestial sphere) in right ascension and declination"""
        return self.data[["ra_tangent", "dec_tangent"]]

    @property
    def radec_bounds(self) -> tuple[float, float, float, float]:
        """bounds in right ascension and declination in order [xmin, ymin, xmax, ymax]"""
        return self.data[["ra_min", "dec_min", "ra_max", "dec_max"]]

    @property
    def orientation(self) -> float:
        return self.data["orientat"]

    @property
    def xy_tangent(self) -> tuple[float, float]:
        """projection origin (tangent point with the celestial sphere) in pixel coordinates"""
        return self.data[["x_tangent", "y_tangent"]]

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels across"""
        return self.data[["nx", "ny"]]

    @cached_property
    def skycell_indices(self) -> NDArray[int]:
        """indices of sky cells in the sky map within this region"""
        return np.arange(self.data["skycell_start"], self.data["skycell_end"])

    @property
    def skycells(self) -> np.void:
        """subset array of sky cells from the sky map within this region"""
        return SKYMAP.skycells[self.skycell_indices]

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
                (self.data["ra_min"], self.data["dec_min"]),
                (self.data["ra_max"], self.data["dec_min"]),
                (self.data["ra_max"], self.data["dec_max"]),
                (self.data["ra_min"], self.data["dec_max"]),
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
        return sga.length(self.vectorpoint_corners[0], self.vectorpoint_corners[2])

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
        return SKYMAP.pixel_scale

    def __eq__(self, other) -> bool:
        if not isinstance(other, ProjectionRegion):
            return NotImplemented

        return self.data == other.data


def wcsinfo_to_wcs(
    wcsinfo: dict | stnode.Wcsinfo,
    bounding_box: None | tuple[tuple[float, float], tuple[float, float]] = None,
    name: str = "wcsinfo",
) -> WCS:
    """Create a WCS from the L3 wcsinfo meta

    Parameters
    ----------
    wcsinfo : dict or MosaicModel.meta.wcsinfo
        The L3 wcsinfo to create a WCS from.

    bounding_box : None or 4-tuple
        The bounding box in detector/pixel space. Form of input is:
        (x_left, x_right, y_bottom, y_top)

    name : str
        Value of the `name` attribute of the WCS object.

    Returns
    -------
    wcs : WCS
        The WCS object created.
    """
    pixelshift = models.Shift(-wcsinfo["x_ref"], name="crpix1") & models.Shift(
        -wcsinfo["y_ref"], name="crpix2"
    )
    pixelscale = models.Scale(wcsinfo["pixel_scale"], name="cdelt1") & models.Scale(
        wcsinfo["pixel_scale"], name="cdelt2"
    )
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        wcsinfo["ra_ref"], wcsinfo["dec_ref"], 180.0
    )

    matrix = wcsinfo.get("rotation_matrix", None)
    if matrix is not None:
        matrix = np.array(matrix)
    else:
        matrix = np.reshape(
            wcs_util.calc_rotation_matrix(
                np.deg2rad(wcsinfo.get("orientat", 0.0)), v3i_yangle=0.0, vparity=1
            ),
            (2, 2),
        )
    rotation = models.AffineTransformation2D(matrix, name="pc_rotation_matrix")
    det2sky = (
        pixelshift | rotation | pixelscale | tangent_projection | celestial_rotation
    )
    det2sky.name = "linear_transform"

    detector_frame = coordinate_frames.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = coordinate_frames.CelestialFrame(
        reference_frame=coordinates.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = WCS([(detector_frame, det2sky), (sky_frame, None)], name=name)

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj


def skycell_to_wcs(skycell_record: dict) -> WCS:
    """From a skycell record, generate a WCS

    Parameters
    ----------
    skycell_record : dict
        A skycell record, or row, from the skycell patches table.

    Returns
    -------
    wcsobj : WCS
        The WCS object from the skycell record.
    """
    wcsinfo = dict()

    # The scale is given in arcseconds per pixel. Convert to degrees.
    wcsinfo["pixel_scale"] = float(skycell_record["pixel_scale"]) / 3600.0

    # Remaining components of the wcsinfo block
    wcsinfo["ra_ref"] = float(skycell_record["ra_projection_center"])
    wcsinfo["dec_ref"] = float(skycell_record["dec_projection_center"])
    wcsinfo["x_ref"] = float(skycell_record["x0_projection"])
    wcsinfo["y_ref"] = float(skycell_record["y0_projection"])
    wcsinfo["orientat"] = float(skycell_record["orientat_projection_center"])

    # Bounding box of the skycell. Note that the center of the pixels are at (0.5, 0.5)
    bounding_box = (
        (-0.5, -0.5 + skycell_record["nx"]),
        (-0.5, -0.5 + skycell_record["ny"]),
    )

    wcsobj = wcsinfo_to_wcs(wcsinfo, bounding_box=bounding_box)

    wcsobj.array_shape = tuple(
        int(axs[1] - axs[0] + 0.5)
        for axs in wcsobj.bounding_box.bounding_box(order="C")
    )
    return wcsobj


def to_skycell_wcs(library: ModelLibrary) -> WCS | None:
    """If available read the skycell WCS from the input library association.

    If the association information contains a "skycell_wcs_info" entry that
    is not "none" it will be interpreted as a skycell wcs. If not, the
    association "target" name will be checked. If it matches a skycell
    name the skycell table will be loaded and a WCS constructed based on the name.
    If neither condition is met None will be returned.

    Parameters
    ----------
    library : ModelLibrary
        ModelLibrary instance containing association information.

    Returns
    -------
    wcsobj : WCS or None
        The WCS object from the skycell record or None if
        none was found.
    """

    try:
        skycell = SkyCell.from_asn(library.asn)

        log.info(f"Skycell record: {skycell.data}")

        # extract the wcs info from the record for skycell_to_wcs
        log.info(
            f"Creating skycell image at ra: {skycell.radec_center[0]}  dec {skycell.radec_center[1]}",
        )

        return skycell.wcs
    except ValueError:
        return None


class SkyMap:
    """
    Abstract representation of the sky map, comprising of 4058 tesellated non-rectangular "sky tiles", each with an area of 10 square degrees.

    For each sky tile, a corresponding gnomonic projection defining a rectangular "projection region" on to a uniform pixel grid entirely covering the sky tile. By necessity, projection regions will overlap other projection regions somewhat. The pixel scale for all projection regions is identical.

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
    def data(self) -> AsdfFile:
        """ASDF representation of sky map"""
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

    @property
    def skycells(self) -> np.void:
        """array of sky cells"""
        return self.data.skycells

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
                        self.skycells["ra_center"], self.skycells["dec_center"]
                    ),
                    axis=1,
                )
            )
        )

    @property
    def projection_regions(self) -> np.void:
        """array of projection regions (one per sky tile)"""
        return self.data.projection_regions

    @cached_property
    def projection_regions_kdtree(self) -> KDTree:
        """k-d tree of all projection regions in the skymap, using normalized center vectorpoints in 3D space"""
        return KDTree(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(
                        self.projection_regions["ra_tangent"],
                        self.projection_regions["dec_tangent"],
                    ),
                    axis=1,
                )
            )
        )

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        # The scale is given in arcseconds per pixel. Convert to degrees.
        return self.data.meta["pixel_scale"] / 3600.0

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels per sky cell"""
        return self.data.meta.nxy_skycell, self.data.meta.nxy_skycell

    def __getitem__(self, index: int) -> SkyCell:
        """`SkyCell` at the given index in the sky cells array"""
        return SkyCell(index)


SKYMAP = SkyMap(path=os.environ.get("SKYMAP_PATH", None))
