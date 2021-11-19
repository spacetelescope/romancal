import numpy as np
from astropy.modeling.models import RotationSequence3D, Scale
from gwcs.geometry import SphericalToCartesian, CartesianToSpherical


def v23tosky(input_model, wrap_v2_at=180, wrap_lon_at=360):
    """ Create the transform from telescope to sky.

    The transform is defined with a reference point in a Frame
    associated tih the telescope (V2, V3) in arcsec, the corresponding
    reference poiont on sky (RA_REF, DEC_REF) in deg, and the position angle
    at the center of the aperture, ROLL_REF in deg.

    Parameters
    ----------
    input_model : `roman_datamodels.WfiImage`
        Roman imaging exposure data model.
    wrap_v2_at : float
        At what angle to wrap V2. [deg]
    wrap_lon_at : float
        At what angle to wrap logitude. [deg]

    Returns
    -------
    model : `astropy.modeling.Model`
        The transform from V2,V3 to sky.
    """
    v2_ref = input_model.meta.wcsinfo.v2_ref / 3600
    v3_ref = input_model.meta.wcsinfo.v3_ref / 3600
    roll_ref = input_model.meta.wcsinfo.roll_ref
    ra_ref = input_model.meta.wcsinfo.ra_ref
    dec_ref = input_model.meta.wcsinfo.dec_ref

    angles = np.array([v2_ref, -v3_ref, roll_ref, dec_ref, -ra_ref])
    axes = "zyxyz"
    rot = RotationSequence3D(angles, axes_order=axes)

    # The sky rotation expects values in deg.
    # This should be removed when models work with quantities.
    model = ((Scale(1 / 3600) & Scale(1 / 3600)) | SphericalToCartesian(wrap_lon_at=wrap_v2_at)
         | rot | CartesianToSpherical(wrap_lon_at=wrap_lon_at))
    model.name = 'v23tosky'
    return model
