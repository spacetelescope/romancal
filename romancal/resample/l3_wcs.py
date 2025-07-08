import logging

import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models
from gwcs import WCS, coordinate_frames
from stcal.alignment.util import (
    compute_s_region_keyword,
    compute_scale,
)

from romancal.assign_wcs.utils import create_footprint

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def assign_l3_wcs(model, wcs):
    """Assign a wcs to a level 3 model and update the wcsinfo

    Parameters
    ----------
    model : `DataModel`
        The model whose meta is to be updated.

    wcs : `GWCS`
        GWCS info to transfer into the `meta.wcsinfo` block

    Notes
    -----
    Some models/parameters in the GWCS object have explicit names, such as
    'crpix1'. However, some do not and hence have to be accessed explicitly
    by indexing. This is fragile and will be a source of issues.
    """
    model.meta.wcs = wcs
    l3_wcsinfo = model.meta.wcsinfo
    transform = wcs.forward_transform

    l3_wcsinfo.projection = "TAN"
    l3_wcsinfo.image_shape = model.shape   # Should this be [::-1], pixel_shape is in cartesian coordinates

    # Fill out image-local information
    pixel_center = [(v - 1) / 2.0 for v in model.shape[::-1]]
    world_center = wcs(*pixel_center)
    l3_wcsinfo.ra = world_center[0]
    l3_wcsinfo.dec = world_center[1]
    l3_wcsinfo.pixel_scale_ref = compute_scale(wcs, world_center)
    l3_wcsinfo.orientation = calc_pa(wcs, *world_center)

    footprint = create_footprint(wcs, model.shape, center=False)
    l3_wcsinfo.s_region = compute_s_region_keyword(footprint)

    try:
        l3_wcsinfo.x_ref = transform.crpix[0]
        l3_wcsinfo.y_ref = transform.crpix[1]
    except AttributeError:
        log.warning(
            "WCS has no clear reference pixel defined by crpix1/crpix2. Assuming reference pixel is center."
        )
        l3_wcsinfo.x_ref = pixel_center[0]
        l3_wcsinfo.y_ref = pixel_center[1]

    world_ref = wcs(l3_wcsinfo.x_ref, l3_wcsinfo.y_ref, with_bounding_box=False)
    l3_wcsinfo.ra_ref = world_ref[0]
    l3_wcsinfo.dec_ref = world_ref[1]

    try:
        cdelt1 = transform.cdelt[0]
        cdelt2 = transform.cdelt[1]
        l3_wcsinfo.pixel_scale = (cdelt1 + cdelt2) / 2.0
    except AttributeError:
        l3_wcsinfo.pixel_scale = compute_scale(wcs, world_ref)

    l3_wcsinfo.orientation_ref = calc_pa(wcs, *world_ref)

    try:
        l3_wcsinfo.rotation_matrix = transform.pc.value.tolist()
    except AttributeError:
        log.warning(
            "WCS has no clear rotation matrix defined by pc_rotation_matrix. Calculating one."
        )
        rotation_matrix = utils.calc_rotation_matrix(l3_wcsinfo.orientat, 0.0)
        l3_wcsinfo.rotation_matrix = utils.list_1d_to_2d(rotation_matrix, 2)


def calc_pa(wcs, ra, dec):
    """Calculate position angle at given ra,dec

    Parameters
    ----------
    wcs : GWCS
        The wcs in consideration.

    ra, dec : float, float
        The ra/dec in degrees.

    Returns
    -------
    position_angle : float
        The position angle in degrees.

    """
    delta_pix = [v for v in wcs.invert(ra, dec, with_bounding_box=False)]
    delta_pix[1] += 100
    delta_coord = SkyCoord(
        *wcs(*delta_pix, with_bounding_box=False), frame="icrs", unit="deg"
    )
    coord = SkyCoord(ra, dec, frame="icrs", unit="deg")

    return coord.position_angle(delta_coord).degree


def l3wcsinfo_to_wcs(wcsinfo, bounding_box=None):
    pixelshift = models.Shift(
        -wcsinfo["x_ref"],
        name="crpix1",
    ) & models.Shift(-wcsinfo["y_ref"], name="crpix2")
    pixelscale = models.Scale(wcsinfo["pixel_scale_ref"], name="scale") & models.Scale(
        wcsinfo["pixel_scale_ref"], name="scale"
    )
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        wcsinfo["ra_ref"], wcsinfo["dec_ref"], 180.0
    )

    rotation = models.AffineTransformation2D(
        np.array(wcsinfo["rotation_matrix"]), name="pc_rotation_matrix"
    )
    det2sky = (
        pixelshift | rotation | pixelscale | tangent_projection | celestial_rotation
    )
    det2sky.name = "detector_to_sky"

    detector_frame = coordinate_frames.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = coordinate_frames.CelestialFrame(
        reference_frame=coordinates.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = WCS(
        [(detector_frame, det2sky), (sky_frame, None)], name=wcsinfo.get("name", None)
    )

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj
