import logging
import os
import re
from datetime import datetime
from functools import cached_property
from pathlib import Path

import asdf
import crds
import gwcs
import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from asdf import AsdfFile
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from astropy.wcs import WCS
from numpy.typing import NDArray
from roman_datamodels import stnode
from scipy.spatial import KDTree
from stcal.alignment import util as wcs_util

from romancal.datamodels.library import ModelLibrary

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


SKYCELL_AREA = (4.6 / 3600.0) ** 2
SKYTILE_AREA = (0.5) ** 2


def image_coords_to_vec(
    radec_coords: list[tuple[float, float]]
    | tuple[list[float], list[float]]
    | NDArray[float],
) -> NDArray[float]:
    """
    This routine can handle the corners in both organizations, whether
    a sequence of ra, dec pairs or a list of ra coordinates, followed
    by a list of dec coordinates. So long as either form can be converted
    to a numpy array it is ok. The spherical geometry routines expect
    the latter organization, and the returned value will be of that
    organization.
    """

    radec_coords = np.array(radec_coords)
    if radec_coords.shape[0] == 2:
        radec_coords = radec_coords.transpose()

    # Convert all celestial coordinates to cartesian coordinates.
    return np.stack(
        sgv.lonlat_to_vector(radec_coords[:, 0], radec_coords[:, 1]), axis=1
    )


class SkyCell:
    """
    square subpartition of a projection region, approx 4.6' x 4.6' and 5000 pixels across
    """

    __index: int | None
    __data: np.void

    def __init__(self, index: int | None):
        """
        :param index: index in the global sky map
        """
        self.__index = index
        if index is not None:
            self.__data = SKYMAP.skycells[index]

    @classmethod
    def from_name(cls, name: str) -> "SkyCell":
        """
        retrieve a sky cell from the sky map by its name


        :param name: name of a sky cell, for instance `315p86x50y75`
        """
        if not re.match(r"\d{3}\w\d{2}x\d{2}y\d{2}", name):
            raise ValueError(f"invalid skycell name {name}")
        return SkyCell(np.where(SKYMAP.skycells["name"] == name)[0][0])

    @classmethod
    def from_data(cls, data: np.void) -> "SkyCell":
        """
        build an index-less sky cell instance from a data array
        :param data: array with sky cell parameters (see schema)
        """
        instance = cls(index=None)
        instance.__data = data
        return instance

    @classmethod
    def from_center_and_coordinates(
        cls,
        skytile_center: tuple[float, float],
        skycell_coordinates: tuple[int, int],
    ) -> "SkyCell":
        """
        retrieve a sky cell from the sky map using its location within its containins sky tile

        :param skytile_center: center coordinates of its containing sky tile in right ascension and declination
        :param skycell_coordinates: XY location of the sky cell within its sky tile, in units of ordinal sky cells from the center
        """
        return cls.from_name(
            f"r{round(skytile_center[0]):03}d{'p' if skytile_center[1] >= 0 else 'm'}{round(skytile_center[1]):02}x{'p' if skycell_coordinates[1] >= 0 else 'm'}{skycell_coordinates[1]:02}y{'p' if skycell_coordinates[1] >= 0 else 'm'}{skycell_coordinates[1]:02}"
        )

    @classmethod
    def from_modellibrary(cls, library: ModelLibrary) -> "SkyCell":
        """
        retrieve a sky cell from WCS info or a target specified in a model library

        :param library: ModelLibrary instance
        """
        if (
            "skycell_wcs_info" in library.asn
            and library.asn["skycell_wcs_info"] != "none"
        ):
            skycell_name = library.asn["skycell_wcs_info"]["name"]
        elif "target" in library.asn:
            # check to see if the product name contains a skycell name & if true get the skycell record
            skycell_name = library.asn["target"]
        else:
            raise ValueError(
                "cannot extract skycell information from modellibrary association with neither WCS nor target info"
            )

        return SkyCell.from_name(skycell_name)

    @property
    def index(self) -> int | None:
        """index of this sky cell in the sky map"""
        return self.__index

    @property
    def data(self) -> np.void:
        """properties from the sky map"""
        return self.__data

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
        """orientation from zenith"""
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

    @property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this sky cell"""
        # convert all radec points to vectors
        corner_vectorpoints = image_coords_to_vec(self.radec_corners)
        center_vectorpoint = image_coords_to_vec([self.radec_center])

        # construct polygon from corner points and center point
        return sgp.SingleSphericalPolygon(
            points=sgv.normalize_vector(corner_vectorpoints),
            inside=sgv.normalize_vector(center_vectorpoint),
        )

    @cached_property
    def projregion(self) -> "ProjectionRegion":
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
        return {
            "ra_ref": self.projregion.data["ra_tangent"],
            "dec_ref": self.projregion.data["dec_tangent"],
            "x_ref": self.projregion.data["x_tangent"],
            "y_ref": self.projregion.data["y_tangent"],
            "rotation_matrix": None,
            "orientat": self.data["orientat"],
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
        """FITS WCS object representing this sky cell"""
        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = [self.pixel_scale, self.pixel_scale]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crval = list(self.projregion.radec_tangent)
        wcs.wcs.crpix = list(self.xy_tangent)
        wcs.wcs.crota = [0, self.projregion.orientation]  # CROTA2 is the rotation angle
        wcs.array_shape = list(self.pixel_shape)
        return wcs

    def __eq__(self, other) -> bool:
        if not isinstance(other, SkyCell):
            return NotImplemented

        return self.data == other.data


class ProjectionRegion:
    """projection region in the sky map, corresponding to a single sky tile"""

    __index: int | None
    __data: np.void

    def __init__(self, index: int | None):
        """
        :param index: index of the projection region in the sky map array
        """
        self.__index = index
        if index is not None:
            self.__data = SKYMAP.projregions[index]

    @classmethod
    def from_data(cls, data: np.void) -> "ProjectionRegion":
        """
        build an index-less projection region instance from a data array
        :param data: array with projection region parameters (see schema)
        """
        instance = cls(index=None)
        instance.__data = data
        return instance

    @classmethod
    def from_skycell_index(cls, index: int) -> "ProjectionRegion":
        """
        :param index: index of the sky cell
        :returns: projection region corresponding to the given sky cell
        """
        for projregion in SKYMAP.projregions:
            if (
                int(projregion["skycell_start"])
                < index
                < int(projregion["skycell_end"])
            ):
                return cls(projregion[0])
        else:
            raise KeyError(
                f"sky cell index {index} not found in any projection regions"
            )

    @property
    def index(self) -> int | None:
        """index in the sky map"""
        return self.__index

    @property
    def data(self) -> np.void:
        """properties from the sky map"""
        return self.__data

    @property
    def radec_tangent(self) -> tuple[float, float]:
        """projection origin (tangent point with the celestial sphere) in right ascension and declination"""
        return self.data["ra_tangent"], self.data["dec_tangent"]

    @property
    def radec_bounds(self) -> tuple[float, float, float, float]:
        """bounds in right ascension and declination in order [xmin, ymin, xmax, ymax]"""
        return (
            self.data["ra_min"],
            self.data["dec_min"],
            self.data["ra_max"],
            self.data["dec_max"],
        )

    @property
    def orientation(self) -> float:
        """orientation from zenith"""
        return self.data["orientat"]

    @property
    def xy_tangent(self) -> tuple[float, float]:
        """projection origin (tangent point with the celestial sphere) in pixel coordinates"""
        return self.data["x_tangent"], self.data["y_tangent"]

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels across"""
        return self.data["nx"], self.data["ny"]

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

    @property
    def polygon(self) -> sgp.SingleSphericalPolygon:
        """spherical polygon representing this region"""
        # convert all radec points to vectors
        corner_vectorpoints = image_coords_to_vec(self.radec_corners)
        center_vectorpoint = np.mean(corner_vectorpoints, axis=0)

        # construct polygon from corner points and center point
        return sgp.SingleSphericalPolygon(
            points=sgv.normalize_vector(corner_vectorpoints),
            inside=sgv.normalize_vector(center_vectorpoint),
        )

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        return SKYMAP.pixel_scale

    @cached_property
    def wcs(self) -> WCS:
        """FITS WCS object representing this region"""
        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = [self.pixel_scale, self.pixel_scale]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crval = list(self.radec_tangent)
        wcs.wcs.crpix = list(self.xy_tangent)
        wcs.wcs.crota = [0, self.orientation]  # CROTA2 is the rotation angle
        wcs.array_shape = list(self.pixel_shape)
        return wcs

    def __eq__(self, other) -> bool:
        if not isinstance(other, ProjectionRegion):
            return NotImplemented

        return self.data == other.data


def wcsinfo_to_gwcs(
    wcsinfo: dict | stnode.Wcsinfo,
    bounding_box: None | tuple[tuple[float, float], tuple[float, float]] = None,
    name: str = "wcsinfo",
) -> gwcs.WCS:
    """Create a GWCS from the L3 wcsinfo meta

    Parameters
    ----------
    wcsinfo : dict or MosaicModel.meta.wcsinfo
        The L3 wcsinfo to create a GWCS from.

    bounding_box : None or 4-tuple
        The bounding box in detector/pixel space. Form of input is:
        (x_left, x_right, y_bottom, y_top)

    name : str
        Value of the `name` attribute of the GWCS object.

    Returns
    -------
    wcs : wcs.GWCS
        The GWCS object created.
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
    if matrix:
        matrix = np.array(matrix)
    else:
        orientat = wcsinfo.get("orientat", 0.0)
        matrix = wcs_util.calc_rotation_matrix(
            np.deg2rad(orientat), v3i_yangle=0.0, vparity=1
        )
        matrix = np.reshape(matrix, (2, 2))
    rotation = models.AffineTransformation2D(matrix, name="pc_rotation_matrix")
    det2sky = (
        pixelshift | rotation | pixelscale | tangent_projection | celestial_rotation
    )

    detector_frame = gwcs.coordinate_frames.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = gwcs.coordinate_frames.CelestialFrame(
        reference_frame=coordinates.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = gwcs.WCS([(detector_frame, det2sky), (sky_frame, None)], name=name)

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj


def skycell_to_gwcs(skycell_record: dict) -> gwcs.WCS:
    """From a skycell record, generate a GWCS

    Parameters
    ----------
    skycell_record : dict
        A skycell record, or row, from the skycell patches table.

    Returns
    -------
    wcsobj : wcs.GWCS
        The GWCS object from the skycell record.
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

    wcsobj = wcsinfo_to_gwcs(wcsinfo, bounding_box=bounding_box)

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
    wcsobj : wcs.GWCS or None
        The GWCS object from the skycell record or None if
        none was found.
    """

    try:
        skycell = SkyCell.from_modellibrary(library)

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
    Abstract representation of the sky map, comprising of
    - 4096 tesellated non-rectangular tiles ("sky tiles"), each with a corresponding rectangular region of gnomonic projection ("projection region")
    - ~2000 square subpartitions ("sky cells", ~8 million in total) for each projection region, each 4.6' or 5000px across
    """

    __path: None | Path

    def __init__(self, path: None | Path | str = None):
        """
        :param path: load sky map from a custom ASDF file (defaults to `skycells` ref on CRDS)
        """
        if path is not None and not isinstance(path, Path):
            path = Path(path)
        self.__path = path

    @property
    def path(self) -> None | Path:
        """location of sky map reference file on filesystem"""
        return self.__path

    @cached_property
    def data(self) -> AsdfFile:
        """ASDF representation of sky map"""
        if self.__path is None:
            rmap = crds.getreferences(
                {
                    "roman.meta.instrument.name": "WFI",
                    "roman.meta.exposure.start_time": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
                },
                reftypes=["skycells"],
                observatory="roman",
            )
            self.__path = Path(rmap["skycells"])
        with asdf.open(self.__path, memmap=True) as file:
            output = file.copy()
        return output

    @property
    def skycells(self) -> np.void:
        """array of sky cells"""
        return self.data["roman"]["skycells"]

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
    def projregions(self) -> np.void:
        """array of projection regions (one per sky tile)"""
        return self.data["roman"]["projection_regions"]

    @cached_property
    def projregions_kdtree(self) -> KDTree:
        """k-d tree of all projection regions in the skymap, using normalized center vectorpoints in 3D space"""
        return KDTree(
            sgv.normalize_vector(
                np.stack(
                    sgv.lonlat_to_vector(
                        self.projregions["ra_tangent"], self.projregions["dec_tangent"]
                    ),
                    axis=1,
                )
            )
        )

    @property
    def pixel_scale(self) -> float:
        """degrees per pixel"""
        # The scale is given in arcseconds per pixel. Convert to degrees.
        return self.data["roman"]["meta"]["pixel_scale"] / 3600.0

    @property
    def pixel_shape(self) -> tuple[int, int]:
        """number of pixels per sky cell"""
        return self.data["roman"]["meta"]["nxy_skycell"], self.data["roman"]["meta"][
            "nxy_skycell"
        ]

    def __getitem__(self, index: int) -> SkyCell:
        """`SkyCell` at the given index in the sky cells array"""
        return SkyCell(self.data["roman"]["skycells"][index])


SKYMAP = SkyMap(path=os.environ.get("SKYMAP_PATH", None))
