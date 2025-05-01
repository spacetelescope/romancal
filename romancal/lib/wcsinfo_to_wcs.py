import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames
from roman_datamodels import stnode
from stcal.alignment import util as wcs_util

__all__ = ["wcsinfo_to_wcs"]


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
        (x_left, x_right, y_bottom, y_top)

    Returns
    -------
    wcs : WCS
        The WCS object created.
    """
    pixelshift = models.Shift(-wcsinfo["x0_projection"], name="offset") & models.Shift(
        -wcsinfo["y0_projection"], name="offset"
    )
    pixelscale = models.Scale(wcsinfo["pixel_scale"], name="scale") & models.Scale(
        wcsinfo["pixel_scale"], name="scale"
    )
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        wcsinfo["ra_projection_center"], wcsinfo["dec_projection_center"], 180.0
    )

    matrix = wcsinfo.get("rotation_matrix", None)
    if matrix is not None:
        matrix = np.array(matrix)
    else:
        matrix = np.reshape(
            wcs_util.calc_rotation_matrix(
                np.deg2rad(wcsinfo.get("orientat_projection_center", 0.0)),
                v3i_yangle=0.0,
                vparity=1,
            ),
            (2, 2),
        )
    rotation = models.AffineTransformation2D(matrix, name="rotation")
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
    wcsobj = WCS(
        [(detector_frame, det2sky), (sky_frame, None)], name=wcsinfo.get("name", None)
    )

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj
