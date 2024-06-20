"""
This module determines which patch files overlap with the given image.

Currently this assumes that the sky projected borders of all images are straight.
"""

import logging
import os
import os.path

import asdf
import gwcs.wcs as wcs
import numpy as np
import spherical_geometry.polygon as sgp
import spherical_geometry.vector as sgv
from spherical_geometry.vector import normalize_vector

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

RAD_TO_ARCSEC = 180.0 / np.pi * 3600.0

PATCH_TABLE = None
CANDIDATE_RADIUS = 0.5  # Radius of circle for which patch centers must lie within
# to be tested for intersection (in degrees)


def load_patch_table(tablepath=None):
    """
    Load the patch table. If no tablepath is supplied the path is obtained
    from the environmental variable PATCH_TABLE_PATH
    """
    global PATCH_TABLE
    if tablepath is None:
        try:
            tablepath = os.environ["PATCH_TABLE_PATH"]
        except KeyError:
            log.error("PATCH_TABLE_PATH environmental variable not found")
            return
    try:
        with asdf.open(tablepath) as af:
            PATCH_TABLE = af.tree["patches"].copy()
    except FileNotFoundError:
        raise FileNotFoundError("Specified patch table file path not found")


def image_coords_to_vec(image_corners):
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


def find_patch_matches(image_corners, image_shape=None):
    """Find patches that the image overlaps with

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
    A sequence of the indices of all patches that overlap the supplied image
    (in the referenced patch table). The indices may be used to obtain all
    necessary information about the patches.
    """

    if PATCH_TABLE is None:
        raise RuntimeError("No patch table has been loaded")
    if isinstance(image_corners, wcs.WCS):
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
    ptab = PATCH_TABLE
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
    # Compute difference vector between image center and patch centers
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


def get_cartesian_corners(patch):
    """
    Cons] truct a the vertex coordinates for a patch from a patch definition suitable
    for plotting the defined region (ra coordinates, dec coordinates). It returns
    coordinates in the cartesian system so that it avoids the distortions in plots
    due to areas approaching the poles.
    """
    dec_corners = np.array(
        (
            patch["dec_corn1"],
            patch["dec_corn2"],
            patch["dec_corn3"],
            patch["dec_corn4"],
            patch["dec_corn1"],
        )
    )
    ra_corners = np.array(
        (
            patch["ra_corn1"],
            patch["ra_corn2"],
            patch["ra_corn3"],
            patch["ra_corn4"],
            patch["ra_corn1"],
        )
    )
    vec_corners = sgv.lonlat_to_vector(ra_corners, dec_corners)
    return vec_corners


def find_closest_tangent_point(patches, image_corners):
    """
    Out of all listed patches, find the closest tangent point to the center
    coordinate of the image.
    """
    # To deal with the corner case, it is necessary to use spherical_geometry
    # to average the the corners properly.
    # Convert image corners to unit vectors.
    vec_im_corners = image_coords_to_vec(image_corners)
    im_center = np.array(normalize_vector(vec_im_corners.mean(axis=1)))
    tangent_point_set = set()
    patch_tangent_points = [
        (patch["ra_projection_center"], patch["dec_projection_center"])
        for patch in patches
    ]
    for tangent_point in patch_tangent_points:
        tangent_point_set.add(tangent_point)
    unique_tangent_points = list(tangent_point_set)
    # Compute distance for each tangent point from im_center
    dist = [
        ((im_center - np.array(sgv.lonlat_to_vector(*tangent_point))) ** 2).sum()
        for tangent_point in unique_tangent_points
    ]
    sorted_dist_indices = sorted(zip(dist, range(len(dist))))
    sorted_tangent_points = [
        unique_tangent_points[sorted_dist[1]] for sorted_dist in sorted_dist_indices
    ]
    closest_tangent_point = np.array(sgv.lonlat_to_vector(*sorted_tangent_points[0]))
    # Now associate index of sorted_tangent_points with that of all patches
    patch_tp_id = []
    for patch_tp in patch_tangent_points:
        for i, tp in enumerate(sorted_tangent_points):
            if tp == patch_tp:
                patch_tp_id.append(i)
    return closest_tangent_point, patch_tp_id


def veccoords_to_tangent_plane(vertices, tangent_point_vec):
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


load_patch_table()
