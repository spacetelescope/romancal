import logging

from astropy.coordinates import SkyCoord
from stcal.alignment.util import (
    #calc_rotation_matrix,
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
    l3_wcsinfo.image_shape = model.shape

    # Fill out image-local information
    pixel_center = [(v - 1) / 2.0 for v in model.shape[::-1]]
    world_center = wcs(*pixel_center)
    l3_wcsinfo.ra = world_center[0]
    l3_wcsinfo.dec = world_center[1]
    l3_wcsinfo.pixel_scale = compute_scale(wcs, world_center)
    l3_wcsinfo.orientation = calc_pa(wcs, *world_center)

    footprint = create_footprint(wcs, model.shape, center=False)
    l3_wcsinfo.s_region = compute_s_region_keyword(footprint)

    # Fill out wcs-general information
    try:
        l3_wcsinfo.x_ref = -transform["crpix1"].offset.value
        l3_wcsinfo.y_ref = -transform["crpix2"].offset.value
    except IndexError:
        log.warning(
            "WCS has no clear reference pixel defined by crpix1/crpix2. Assuming reference pixel is center."
        )
        l3_wcsinfo.x_ref = pixel_center[0]
        l3_wcsinfo.y_ref = pixel_center[1]

    world_ref = wcs(l3_wcsinfo.x_ref, l3_wcsinfo.y_ref, with_bounding_box=False)
    l3_wcsinfo.ra_ref = world_ref[0]
    l3_wcsinfo.dec_ref = world_ref[1]

    try:
        cdelt1 = transform["cdelt1"].factor.value
        cdelt2 = transform["cdelt2"].factor.value
        l3_wcsinfo.pixel_scale_ref = (cdelt1 + cdelt2) / 2.0
    except IndexError:
        l3_wcsinfo.pixel_scale_ref = compute_scale(wcs, world_ref)

    l3_wcsinfo.orientation_ref = calc_pa(wcs, *world_ref)

    #try:
    #    l3_wcsinfo.rotation_matrix = transform[
    #        "pc_rotation_matrix"
    #    ].matrix.value.tolist()
    #except Exception:
    #    log.warning(
    #        "WCS has no clear rotation matrix defined by pc_rotation_matrix. Calculating one."
    #    )
    #    rotation_matrix = calc_rotation_matrix(l3_wcsinfo.orientat, 0.0)
    #    l3_wcsinfo.rotation_matrix = _list_1d_to_2d(rotation_matrix, 2)


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


# def _list_1d_to_2d(l, n):
#     """Convert 1-dimensional list to 2-dimensional
# 
#     Parameters
#     ----------
#     l : list
#         The list to convert.
# 
#     n : int
#        The length of the x dimension, or the length of the inner lists.
# 
#     Returns
#     -------
#     l2d : list of lists
#         The 2D form
#     """
#     return [l[i : i + n] for i in range(0, len(l), n)]
