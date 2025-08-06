import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames, fitswcs
from roman_datamodels import stnode
from stcal.alignment import util as wcs_util

__all__ = ["wcsinfo_to_wcs"]


def wcsinfo_to_wcs(
    wcsinfo: dict | stnode.Wcsinfo,
    bounding_box: None | tuple[tuple[float, float], tuple[float, float]] = None,
) -> WCS:
    """Create a WCS from the skycell wcsinfo meta

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
    crpix = [wcsinfo["x0_projection"], wcsinfo["y0_projection"]]
    crval = [wcsinfo["ra_projection_center"], wcsinfo["dec_projection_center"]]
    cdelt = [wcsinfo.get("pixel_scale", 1), wcsinfo.get("pixel_scale", 1)]

    tangent_projection = models.Pix2Sky_TAN()

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

    det2sky = fitswcs.FITSImagingWCSTransform(
        tangent_projection, crpix=crpix, crval=crval, cdelt=cdelt, pc=matrix
    )

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
