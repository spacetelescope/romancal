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
    """Create a WCS from the L3 wcsinfo meta or a skycell reference file entry

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
    # this has to handle multiple formats:
    # -- skycells ----------------- L3
    # - name                        -
    # - pixel_scale                 - pixel_scale_ref
    # - ra_projection_center        - (ra_ref)
    # - dec_projection_center       - (dec_ref)
    # - x0_projection               - (x_ref)
    # - y0_projection               - (y_ref)
    # - ra_center                   - ra
    # - dec_center                  - dec
    # - nx                          -
    # - ny                          -
    # - orientat                    - orientation_ref
    # - orientat_projection_center  -
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

    # no input has a rotation matrix
    matrix = np.reshape(
        wcs_util.calc_rotation_matrix(
            np.deg2rad(
                wcsinfo.get(
                    "orientat_projection_center", wcsinfo.get("orientation_ref", 0.0)
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

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj
