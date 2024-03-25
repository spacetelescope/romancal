"""
This module determines which tessel files overlap with the given image.

Currently this assumes that the sky projected borders of all images are straight.
"""

import os
import os.path
import numpy as np
import asdf
import spherical_geometry.vector as sgv
import spherical_geometry.polygon as sgp
import gwcs.wcs as wcs
from matplotlib import pyplot as plt

RAD_TO_ARCSEC = 180. / np.pi * 3600.

plt.ion()

TESSEL_TABLE = None

print(os.environ)

def load_tessel_table(tablepath=None):
    """
    Load the tessel table. If no tablepath is supplied the path is obtained
    from the environmental variable TESSEL_TABLE_PATH
    """
    global TESSEL_TABLE
    if tablepath is None:
        try:
            tablepath = os.environ['TESSEL_TABLE_PATH']
        except KeyError:
            raise KeyError("TESSEL_TABLE_PATH environmental variable not found")
    try:
        with asdf.open(tablepath) as af:
            TESSEL_TABLE = af.tree['tessellation'].copy()
    except FileNotFoundError:
        raise FileNotFoundError("Specified tessel table file path not found")


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

def find_tessel_matches(image_corners, image_shape=None):
    """Find tessels that the image overlaps with

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
    A sequence of the indices of all tessels that overlap the supplied image
    (in the referenced tessel table). The indices may be used to obtain all
    necessary information about the tessels.
    """

    if isinstance(image_corners, wcs.WCS):
        iwcs = image_corners
        # Now must find size of correspinding image, with three possible
        # sources of that information.
        if ((not hasattr(iwcs, "bounding_box") or iwcs.bounding_box is None)
             and (not hasattr(iwcs, "pixel_shape") or iwcs.pixel_shape is None)
             and image_shape is None):
            raise ValueError(
                "Use of a wcs object requires at least one of the bounding_box"
                " or pixel_shape attributes be set to the image shape or that the"
                "image_shape argument be set")
        if image_shape is not None:
            pass
        else:
            if hasattr(iwcs, "bounding_box") and iwcs.bounding_box is not None:
                # Presumes that the bounding_box matches image array boundaries
                bbintervals = iwcs.bounding_box.intervals
                # This compensates for the half pixel adjustment in the general code.
                image_shape = (bbintervals[0].upper + 0.5, bbintervals[1].upper + 0.5)
            elif hasattr(iwcs, "pixel_shape") and iwcs.pixel_shape is not None:
                image_shape = iwcs.pixel_shape
        # Compute the image corners ra, dec from the wcs
        (cxm, cxp), (cym, cyp) = ((-0.5, image_shape[1] - 0.5),
                                  (-0.5, image_shape[0] - 0.5))
        image_corners = (iwcs(cxp, cyp), iwcs(cxm, cyp), iwcs(cxm, cym), iwcs(cxp, cym))
    ttab = TESSEL_TABLE
    ra = ttab[:]['ra_center']
    dec = ttab[:]['dec_center']
    # # Convert all celestial coordinates to cartesion coordinates.
    vec_centers = np.array(sgv.lonlat_to_vector(ra, dec)).transpose()
    # # Organize corners into two ra, dec lists
    # imra = np.array([point[0] for point in image_corners])
    # imdec = np.array([point[1] for point in image_corners])
    vec_im_corners = image_coords_to_vec(image_corners)
    # Approximate center of image by averaging corner vectors
    im_center = normalize_vector(vec_im_corners.mean(axis=1))
    # Compute difference vector between image centerand tessel centers
    diff = vec_centers - im_center
    dist = np.sqrt((diff**2).sum(axis=1)) * 180 / np.pi
    match = np.where(dist < 0.5)
    ncandidates = len(match[0])
    # Now see which of these that are close actually overlap the supplied image.
    # (Is it necessary to check that the corners are in a sensible order?)
    mra1 = ttab[match]['ra_corn1']
    mra2 = ttab[match]['ra_corn2']
    mra3 = ttab[match]['ra_corn3']
    mra4 = ttab[match]['ra_corn4']
    mdec1 = ttab[match]['dec_corn1']
    mdec2 = ttab[match]['dec_corn2']
    mdec3 = ttab[match]['dec_corn3']
    mdec4 = ttab[match]['dec_corn4']
    mcenters = vec_centers[match]
    mra = np.vstack([mra1, mra2, mra3, mra4, mra1])
    mdec = np.vstack([mdec1, mdec2, mdec3, mdec4, mdec1])
    points = np.array(sgv.lonlat_to_vector(mra, mdec))
    # Create polygon for supplied image_corners
    pvec_im_corners = np.zeros((5, 3))
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

def get_corners(tessel):
    """
    Construct a the vertex coordinates for a tessel from a tessel definition suitable
    for plotting the defined region (x coordinates, y coordinates).
    """
    t = tessel
    corners = ((t['dec_corn1'],
                t['dec_corn2'],
                t['dec_corn3'],
                t['dec_corn4'],
                t['dec_corn1']),
               (t['ra_corn1'],
                t['ra_corn2'],
                t['ra_corn3'],
                t['ra_corn4'],
                t['ra_corn1']))
    corners = np.array(corners)
    vec_corners = sgv.lonlat_to_vector(corners[0], corners[1])
    return vec_corners

def find_closest_tangent_point(tessels, image_corners):
    """
    Out of all listed tessels, find the closest tangent point to the center
    coordinate of the image.
    """
    # To deal with the corner case, it is necessary to use spherical_geometery
    # to average the the corners properly.
    # Convert image corners to unit vectors.
    vec_im_corners = image_coords_to_vec(image_corners)
    im_center = np.array(normalize_vector(vec_im_corners.mean(axis=1)))
    tangent_point_set = set()
    tessel_tangent_points = [(tess['dec_projection_center'],
                              tess['ra_projection_center'])
        for tess in tessels]
    for tangent_point in tessel_tangent_points:
        tangent_point_set.add(tangent_point)
    unique_tangent_points = list(tangent_point_set)
    # Compute distance for each tangent point from im_center
    dist = [((im_center - np.array(sgv.lonlat_to_vector(*tangent_point)))**2).sum()
           for tangent_point in unique_tangent_points]
    sorted_dist_indices = sorted(zip(dist, range(len(dist))))
    sorted_tangent_points = [unique_tangent_points[sorted_dist[1]]
           for sorted_dist in sorted_dist_indices]
    closest_tangent_point = np.array(sgv.lonlat_to_vector(*sorted_tangent_points[0]))
    # Now associate index of sorted_tangent_points with that of all tessels
    tessel_tp_id = []
    for tess_tp in tessel_tangent_points:
        for i, tp in enumerate(sorted_tangent_points):
            if tp == tess_tp:
                tessel_tp_id.append(i)
    return closest_tangent_point, tessel_tp_id

def normalize_vector(vec):
    """
    Normalize a 3d vector to have length 1. Only works on 1d arrays.
    """
    return vec/np.sqrt((vec**2).sum())


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
    y_axis = normalize_vector(np.array([0, 0, 1])
        - np.array(tangent_point_vec) * np.dot(np.array([0, 0, 1]),
        np.array(tangent_point_vec)))
    avertices = np.vstack(vertices)
    x_coords = np.dot(x_axis, avertices) * RAD_TO_ARCSEC
    y_coords = np.dot(y_axis, avertices) * RAD_TO_ARCSEC
    return x_coords, y_coords

def plot_field(corners, id='', fill=None, color=None):
    plt.fill(corners[0], corners[1], color=fill, edgecolor=color)

def plot_tessel(corners, id='', color=None):
    plt.plot(corners[0], corners[1], color=color)
    if id:
        center = (corners[0][:-1].mean(), corners[1][:-1].mean())
        plt.annotate(str(id), center, va='center', ha='center', size=10)

def plot(image_corners, tessels_touched_ids, tessels_candidate_ids):
    """
    This plots a list of tessels
    """
    plt.clf()
    plt.gca().invert_xaxis()
    plt.plot(0, 0, '*', markersize=10)
    tessels_touched = [TESSEL_TABLE[index] for index in tessels_touched_ids]
    tessels_candidate = [TESSEL_TABLE[index] for index in tessels_candidate_ids]
    tangent_point, tess_tp_id_touched = find_closest_tangent_point(
        tessels_touched, image_corners)
    ra, dec = sgv.vector_to_lonlat(*tangent_point)
    dummy, tess_tp_id = find_closest_tangent_point(tessels_candidate, image_corners)
    vec_image_corners = image_coords_to_vec(image_corners)
    tp_image_corners = veccoords_to_tangent_plane(vec_image_corners, tangent_point)
    plot_field(tp_image_corners, fill='lightgrey', color='black')
    for tess, id in zip(tessels_candidate, tessels_candidate_ids):
        plot_tessel(veccoords_to_tangent_plane(get_corners(tess), tangent_point),
            id=id, color='lightgray')
    for tess, id in zip(tessels_touched, tessels_touched_ids):
        plot_tessel(veccoords_to_tangent_plane(get_corners(tess), tangent_point),
            id=id, color='blue')
    plt.xlabel("Offset from nearest tangent point in arcsec")
    plt.ylabel("Offset from nearest tangent point in arcsec")
    plt.title(f"RA: {ra} Dec: {dec} of tangent point in degrees")

load_tessel_table()
