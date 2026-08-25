"""
Helpers for building WCS objects in tests.
"""

from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf

__all__ = ["create_wcs_object"]


def create_wcs_object(fiducial_world, pscale, shape, distortion=None, **kwargs):
    """
    Create a simple tangent-plane WCS object without rotation.

    The fiducial world coordinate lies at pixel ``(0, 0)``.

    Parameters
    ----------
    fiducial_world : tuple
        A pair of values corresponding to the fiducial's world coordinate.
    pscale : tuple
        A pair of values corresponding to the pixel scale in each axis.
    shape : tuple
        A pair of values specifying the dimensions of the WCS object.
    distortion : `~astropy.modeling.Model`, optional
        A two-input, two-output model applied in pixel coordinates
        before the plate scale. If `None`, the WCS has no distortion and
        the pixel area is uniform.

    Returns
    -------
    gwcs.WCS
        A gwcs.WCS object.
    """
    # components of the model
    shift = models.Shift() & models.Shift()
    affine = models.AffineTransformation2D(
        matrix=[[1, 0], [0, 1]], translation=[0, 0], name="pc_rotation_matrix"
    )
    scale = models.Scale(pscale[0]) & models.Scale(pscale[1])
    tan = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        fiducial_world[0],
        fiducial_world[1],
        180,
    )

    # transforms between frames
    # detector -> sky
    det2sky = shift
    if distortion is not None:
        det2sky = det2sky | distortion
    det2sky = det2sky | affine | scale | tan | celestial_rotation
    det2sky.name = "linear_transform"

    # frames
    detector_frame = cf.Frame2D(
        name="detector",
        axes_order=(0, 1),
        axes_names=("x", "y"),
        unit=(u.pix, u.pix),
    )
    sky_frame = cf.CelestialFrame(
        reference_frame=coord.FK5(), name="fk5", unit=(u.deg, u.deg)
    )

    pipeline = [
        (detector_frame, det2sky),
        (sky_frame, None),
    ]

    wcs_obj = WCS(pipeline)

    wcs_obj.bounding_box = kwargs.get(
        "bounding_box",
        (
            (-0.5, shape[-1] - 0.5),
            (-0.5, shape[-2] - 0.5),
        ),
    )

    wcs_obj.pixel_shape = shape[::-1]
    wcs_obj.array_shape = shape

    return wcs_obj
