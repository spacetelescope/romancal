"""
This module determines which projection regions ("sky tiles") overlap with the given image.

Currently this assumes that the sky projected borders of all images are straight.
"""

import logging
import os
import os.path
import re

import asdf
import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames
from numpy.typing import NDArray
from roman_datamodels import stnode
from spherical_geometry.vector import normalize_vector
from stcal.alignment import util as wcs_util

from romancal.datamodels.library import ModelLibrary

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

RAD_TO_ARCSEC = 180.0 / np.pi * 3600.0

PROJREGION_TABLE: dict = None
CANDIDATE_RADIUS = (
    0.5  # Radius of circle for which projection region centers must lie within
)
# to be tested for intersection (in degrees)


def load_projregion_table(tablepath: str | os.PathLike | None = None):
    """
    Load the table of projection regions. If no tablepath is supplied the path is obtained
    from the environmental variable PROJREGION_TABLE_PATH
    """
    global PROJREGION_TABLE
    if tablepath is None:
        try:
            tablepath = os.environ["PROJREGION_TABLE_PATH"]
        except KeyError:
            log.error("PROJREGION_TABLE_PATH environmental variable not found")
            return
    try:
        with asdf.open(tablepath) as af:
            PROJREGION_TABLE = af.tree["roman"]["projection_regions"].copy()
    except FileNotFoundError as err:
        raise FileNotFoundError(
            "Specified projection region table file path not found"
        ) from err


def image_coords_to_vec(
    image_corners: list[tuple[float, float]] | tuple[list[float], list[float]],
) -> NDArray[float]:
    """
    This routine can handle the corners in both organizations, whether
    a sequence of ra, dec pairs or a list of ra coordinates, followed
    by a list of dec coordinates. So long as either form can be converted
    to a numpy array it is ok. The spherical geometry routines expect
    the latter organization, and the returned value will be of that
    organization.
    """
    image_corners = np.array(image_corners)
    if image_corners.shape == (4, 2):
        image_corners = image_corners.transpose()
    # Convert all celestial coordinates to cartesion coordinates.
    vec_im_corners = np.array(sgv.lonlat_to_vector(image_corners[0], image_corners[1]))
    return vec_im_corners


def find_proj_matches(
    image_corners: list[tuple[float, float]] | NDArray[float] | WCS,
    image_shape: tuple[int, int] | None = None,
):
    """Find projection regions that the image overlaps with

    Parameters
    ----------
    image_corners : Either a squence of 4 (ra, dec) pairs, or
        equivalent 2-d numpy array, or
        a GWCS instance. The instance must have either the bounding_box or
        pixel_shape attribute defined, or the following image_shape argument
        must be supplied

    image_shape : image shape to be used if a GWCS instance is supplied
        and does not contain a value for either the bounding_box or
        pixel_shape attributes. Default value is None.

    Returns
    -------
    A sequence of the indices of all projection regions that overlap the supplied image
    (in the referenced projection region table). The indices may be used to obtain all
    necessary information about the projection regions.
    """

    if PROJREGION_TABLE is None:
        load_projregion_table()
    if PROJREGION_TABLE is None:
        raise RuntimeError("No projection region table has been loaded")
    if isinstance(image_corners, WCS):
        iwcs = image_corners
        # Now must find size of correspinding image, with three possible
        # sources of that information.
        if (
            (not hasattr(iwcs, "bounding_box") or iwcs.bounding_box is None)
            and (not hasattr(iwcs, "pixel_shape") or iwcs.pixel_shape is None)
            and image_shape is None
        ):
            raise ValueError(
                "Use of a wcs object requires at least one of the bounding_box"
                " or pixel_shape attributes be set to the image shape or that the"
                "image_shape argument be set"
            )
        if image_shape is not None:
            pass
        else:
            # Both bounding_box and pixel_shape are in x, y order contrary to
            # numpy convention.
            if hasattr(iwcs, "bounding_box") and iwcs.bounding_box is not None:
                # Presumes that the bounding_box matches image array boundaries
                bbintervals = iwcs.bounding_box.intervals
                # This compensates for the half pixel adjustment in the general code.
                image_shape = (bbintervals[1].upper + 0.5, bbintervals[0].upper + 0.5)
            elif hasattr(iwcs, "pixel_shape") and iwcs.pixel_shape is not None:
                image_shape = (iwcs.pixel_shape[1], iwcs.pixel_shape[0])
        # Compute the image corners ra, dec from the wcs
        (cxm, cxp), (cym, cyp) = (
            (-0.5, image_shape[1] - 0.5),
            (-0.5, image_shape[0] - 0.5),
        )
        image_corners = (iwcs(cxp, cyp), iwcs(cxm, cyp), iwcs(cxm, cym), iwcs(cxp, cym))
    ptab = PROJREGION_TABLE
    ra = ptab[:]["ra_center"]
    dec = ptab[:]["dec_center"]
    # # Convert all celestial coordinates to cartesion coordinates.
    vec_centers = np.array(sgv.lonlat_to_vector(ra, dec)).transpose()
    # # Organize corners into two ra, dec lists
    # imra = np.array([point[0] for point in image_corners])
    # imdec = np.array([point[1] for point in image_corners])
    vec_im_corners = image_coords_to_vec(image_corners)
    # Approximate center of image by averaging corner vectors
    im_center = normalize_vector(vec_im_corners.mean(axis=1))
    # Compute difference vector between image center and projection region centers
    diff = vec_centers - im_center
    dist = np.sqrt((diff**2).sum(axis=1)) * 180 / np.pi
    match = np.where(dist < 0.5)
    ncandidates = len(match[0])
    # Now see which of these that are close actually overlap the supplied image.
    # (Is it necessary to check that the corners are in a sensible order?)
    # All the corner coordinates are returned as arrays.
    mra1 = ptab[match]["ra_corn1"]
    mra2 = ptab[match]["ra_corn2"]
    mra3 = ptab[match]["ra_corn3"]
    mra4 = ptab[match]["ra_corn4"]
    mdec1 = ptab[match]["dec_corn1"]
    mdec2 = ptab[match]["dec_corn2"]
    mdec3 = ptab[match]["dec_corn3"]
    mdec4 = ptab[match]["dec_corn4"]
    mcenters = vec_centers[match]
    mra = np.vstack([mra1, mra2, mra3, mra4, mra1])
    mdec = np.vstack([mdec1, mdec2, mdec3, mdec4, mdec1])
    points = np.array(sgv.lonlat_to_vector(mra, mdec))
    # Create polygon for supplied image_corners
    pvec_im_corners = np.empty((5, 3))
    pvec_im_corners[:4] = vec_im_corners.transpose()
    pvec_im_corners[4] = pvec_im_corners[0]
    impoly = sgp.SingleSphericalPolygon(pvec_im_corners, im_center)
    realmatch = []
    for i in range(ncandidates):
        # print(i)
        cellpoints = points[:, :, i].transpose()
        cellcenter = mcenters[i]
        cellpoly = sgp.SingleSphericalPolygon(cellpoints, cellcenter)
        if impoly.intersects_poly(cellpoly):
            # print(f"candidate {i} intersects")
            realmatch.append(i)
    return match[0][realmatch], match[0]


def get_cartesian_corners(
    projection_region: dict[str, NDArray[float]],
) -> tuple[NDArray[float], NDArray[float]]:
    """
    Construct vertex coordinates for a projection region definition suitable
    for plotting the defined region (ra coordinates, dec coordinates). It returns
    coordinates in the cartesian system so that it avoids the distortions in plots
    due to areas approaching the poles.
    """
    dec_corners = np.array(
        (
            projection_region["dec_corn1"],
            projection_region["dec_corn2"],
            projection_region["dec_corn3"],
            projection_region["dec_corn4"],
            projection_region["dec_corn1"],
        )
    )
    ra_corners = np.array(
        (
            projection_region["ra_corn1"],
            projection_region["ra_corn2"],
            projection_region["ra_corn3"],
            projection_region["ra_corn4"],
            projection_region["ra_corn1"],
        )
    )
    vec_corners = sgv.lonlat_to_vector(ra_corners, dec_corners)
    return vec_corners


def find_closest_tangent_point(
    projection_regions: list[dict],
    image_corners: list[tuple[float, float]] | tuple[list[float], list[float]],
):
    """
    Out of all listed projection regions, find the closest tangent point to the center
    coordinate of the image.
    """
    # To deal with the corner case, it is necessary to use spherical_geometry
    # to average the the corners properly.
    # Convert image corners to unit vectors.
    vec_im_corners = image_coords_to_vec(image_corners)
    im_center = np.array(normalize_vector(vec_im_corners.mean(axis=1)))
    tangent_point_set = set()
    proj_tangent_points = [
        (region["ra_projection_center"], region["dec_projection_center"])
        for region in projection_regions
    ]
    for proj_tangent_point in proj_tangent_points:
        tangent_point_set.add(proj_tangent_point)
    unique_tangent_points = list(tangent_point_set)
    # Compute distance for each tangent point from im_center
    dist = [
        ((im_center - np.array(sgv.lonlat_to_vector(*tangent_point))) ** 2).sum()
        for tangent_point in unique_tangent_points
    ]
    sorted_dist_indices = sorted(zip(dist, range(len(dist)), strict=False))
    sorted_tangent_points = [
        unique_tangent_points[sorted_dist[1]] for sorted_dist in sorted_dist_indices
    ]
    closest_tangent_point = np.array(sgv.lonlat_to_vector(*sorted_tangent_points[0]))
    # Now associate index of sorted_tangent_points with that of all projection regions
    projregion_tangentpoint_indices = [
        index
        for proj_tangent_point in proj_tangent_points
        for index, tangent_point in enumerate(sorted_tangent_points)
        if tangent_point == proj_tangent_point
    ]
    return closest_tangent_point, projregion_tangentpoint_indices


def veccoords_to_tangent_plane(
    vertices: list[tuple[float, float, float]],
    tangent_point_vec: tuple[float, float, float],
):
    """
    Convert the spherical geometry vectors to tangent plane coordinates
    in arcseconds. This algorithm is not precise, but should be good
    enough for now (and besides, the goal here is visualizaion, not
    ultra-precision). This also breaks down numerically very near the
    poles.
    """
    # First compute the tangent plane axis vectors.
    x_axis = normalize_vector(np.cross([0, 0, 1], tangent_point_vec))
    y_axis = normalize_vector(
        np.array([0, 0, 1])
        - np.array(tangent_point_vec)
        * np.dot(np.array([0, 0, 1]), np.array(tangent_point_vec))
    )
    avertices = np.vstack(vertices)
    x_coords = np.dot(x_axis, avertices) * RAD_TO_ARCSEC
    y_coords = np.dot(y_axis, avertices) * RAD_TO_ARCSEC
    return x_coords, y_coords


def wcsinfo_to_wcs(
    wcsinfo: dict | stnode.Wcsinfo,
    bounding_box: None | tuple[float, float, float, float] = None,
    name: str = "wcsinfo",
) -> WCS:
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


def skycell_to_wcs(skycell_record):
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
    wcsinfo["rotation_matrix"] = None

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


def to_skycell_wcs(library: ModelLibrary) -> WCS:
    """If available read the skycell WCS from the input library association.

    If the association information contains a "skycell_wcs_info" entry that
    is not "none" it will be interpreted as a skycell wcs. If not, the
    association "target" name will be checked. If it matches a skycell
    name the patch table will be loaded and the skycell wcs will be
    looked up based on the name. If neither condition is met None
    will be returned.

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

    if "skycell_wcs_info" in library.asn and library.asn["skycell_wcs_info"] != "none":
        skycell_record = library.asn["skycell_wcs_info"]
    else:
        if "target" not in library.asn:
            return None
        # check to see if the product name contains a skycell name & if true get the skycell record
        skycell_name = library.asn["target"]

        if not re.match(r"r\d{3}\w{2}\d{2}x\d{2}y\d{2}", skycell_name):
            return None

        if PROJREGION_TABLE is None:
            load_projregion_table()
        skycell_record = PROJREGION_TABLE[
            np.where(PROJREGION_TABLE["name"][:] == skycell_name)[0][0]
        ]
    log.info("Skycell record %s:", skycell_record)

    # extract the wcs info from the record for skycell_to_wcs
    log.info(
        "Creating skycell image at ra: %f  dec %f",
        float(skycell_record["ra_center"]),
        float(skycell_record["dec_center"]),
    )
    return skycell_to_wcs(skycell_record)


def get_projectioncell_wcs(index: np.int64) -> dict | None:
    """Return the projection cell wcs info as a dictionary based on the db index number"""

    # check to see if an index is being passed
    if not isinstance(index, np.int64):
        log.info("Input index needs to be a numpy int64 variable")
        return None

    # check to see if the patch table is loaded if not load it
    if PROJREGION_TABLE is None:
        load_projregion_table()

    projcell_info = dict(
        [
            ("name", PROJREGION_TABLE[index]["name"]),
            ("pixel_scale", float(PROJREGION_TABLE[index]["pixel_scale"])),
            (
                "ra_projection_center",
                float(PROJREGION_TABLE[index]["ra_projection_center"]),
            ),
            (
                "dec_projection_center",
                float(PROJREGION_TABLE[index]["dec_projection_center"]),
            ),
            ("x0_projection", float(PROJREGION_TABLE[index]["x0_projection"])),
            ("y0_projection", float(PROJREGION_TABLE[index]["y0_projection"])),
            ("ra_center", float(PROJREGION_TABLE[index]["ra_center"])),
            ("dec_center", float(PROJREGION_TABLE[index]["dec_center"])),
            ("nx", int(PROJREGION_TABLE[index]["nx"])),
            ("ny", int(PROJREGION_TABLE[index]["ny"])),
            ("orientat", float(PROJREGION_TABLE[index]["orientat"])),
            (
                "orientat_projection_center",
                float(PROJREGION_TABLE[index]["orientat_projection_center"]),
            ),
        ]
    )

    return projcell_info
