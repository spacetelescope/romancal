import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames
from roman_datamodels import stnode
from stcal.alignment import util as wcs_util

__all__ = ["wcsinfo_to_wcs"]


def _skycell_wcsinfo_to_wcs(wcsinfo):
    pixelshift = models.Shift(
        -wcsinfo.get("x0_projection", wcsinfo.get("x_ref", None)),
        name="crpix1",
    ) & models.Shift(
        -wcsinfo.get("y0_projection", wcsinfo.get("y_ref", None)),
        name="crpix2",
    )
    pixelscale = models.Scale(wcsinfo["pixel_scale"], name="scale") & models.Scale(
        wcsinfo["pixel_scale"], name="scale"
    )
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        wcsinfo.get("ra_projection_center", wcsinfo.get("ra_ref", None)),
        wcsinfo.get("dec_projection_center", wcsinfo.get("dec_ref", None)),
        180.0,
    )

    matrix = wcsinfo.get("rotation_matrix", None)
    if matrix is not None:
        matrix = np.array(matrix)
    else:
        matrix = np.reshape(
            wcs_util.calc_rotation_matrix(
                np.deg2rad(
                    wcsinfo.get(
                        "orientat_projection_center", wcsinfo.get("orientat", 0.0)
                    )
                ),
                v3i_yangle=0.0,
                vparity=1,
            ),
            (2, 2),
        )
    rotation = models.AffineTransformation2D(matrix, name="pc_rotation_matrix")
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

    return wcsobj


# mapping from L3 wcsinfo metadata keywords to
# keywords used for the skycell wcsinfo
_L3_TO_SKYCELL_MAPPING = {
    "skycell_name": "name",
    "pixel_scale_ref": "pixel_scale",
    "ra_ref": "ra_projection_center",
    "dec_ref": "dec_projection_center",
    "x_ref": "x0_projection",
    "y_ref": "y0_projection",
    "ra": "ra_center",
    "dec": "dec_center",
    "orietation_ref": "orientat",
}


def _l3_wcsinfo_to_wcs(wcsinfo):
    unmapped_keys = set(wcsinfo.keys())
    mapped_info = {}
    for from_key, to_key in _L3_TO_SKYCELL_MAPPING.items():
        if from_key in wcsinfo:
            unmapped_keys.discard(from_key)
            unmapped_keys.discard(to_key)
            mapped_info[to_key] = wcsinfo[from_key]
    for key in unmapped_keys:
        mapped_info[key] = wcsinfo[key]
    return _skycell_wcsinfo_to_wcs(mapped_info)


def wcsinfo_to_wcs(
    wcsinfo: dict | stnode.Wcsinfo,
    bounding_box: None | tuple[tuple[float, float], tuple[float, float]] = None,
) -> WCS:
    """Create a WCS from the L3 wcsinfo meta

    Parameters
    ----------
    wcsinfo : dict or MosaicModel.meta.wcsinfo
        The L3 wcsinfo to create a WCS from.

    bounding_box : None or 4-tuple
        The bounding box in detector/pixel space. Form of input is:
        ((x_left, x_right), (y_bottom, y_top))

    Returns
    -------
    wcs : WCS
        The WCS object created.
    """
    if "orientat" in wcsinfo or "rotation_matrix" in wcsinfo:
        wcsobj = _skycell_wcsinfo_to_wcs(wcsinfo)
    else:
        wcsobj = _l3_wcsinfo_to_wcs(wcsinfo)

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj
