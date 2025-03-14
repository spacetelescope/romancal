import logging

import numpy as np
from astropy.modeling.models import Identity, RotationSequence3D, Scale, Shift
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def v23tosky(input_model, wrap_v2_at=180, wrap_lon_at=360):
    """Create the transform from telescope to sky.

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
    model = (
        (Scale(1 / 3600) & Scale(1 / 3600))
        | SphericalToCartesian(wrap_lon_at=wrap_v2_at)
        | rot
        | CartesianToSpherical(wrap_lon_at=wrap_lon_at)
    )
    model.name = "v23tosky"
    return model


def dva_corr_model(va_scale, v2_ref, v3_ref):
    """
    Create transformation that accounts for differential velocity aberration
    (scale).

    Parameters
    ----------
    va_scale : float, None
        Ratio of the apparent plate scale to the true plate scale. When
        ``va_scale`` is `None`, it is assumed to be identical to ``1`` and
        an ``astropy.modeling.models.Identity`` model will be returned.

    v2_ref : float, None
        Telescope ``v2`` coordinate of the reference point in ``arcsec``. When
        ``v2_ref`` is `None`, it is assumed to be identical to ``0``.

    v3_ref : float, None
        Telescope ``v3`` coordinate of the reference point in ``arcsec``. When
        ``v3_ref`` is `None`, it is assumed to be identical to ``0``.

    Returns
    -------
    va_corr : astropy.modeling.CompoundModel, astropy.modeling.models.Identity
        A 2D compound model that corrects DVA. If ``va_scale`` is `None` or 1
        then `astropy.modeling.models.Identity` will be returned.

    """
    if va_scale is None or va_scale == 1:
        return Identity(2)

    if va_scale <= 0:
        log.warning("Given velocity aberration scale %s", va_scale)
        log.warning(
            "Velocity aberration scale must be a positive number. Setting to 1.0"
        )
        va_scale = 1.0

    va_corr = Scale(va_scale, name="dva_scale_v2") & Scale(
        va_scale, name="dva_scale_v3"
    )

    if v2_ref is None:
        v2_ref = 0

    if v3_ref is None:
        v3_ref = 0

    if v2_ref == 0 and v3_ref == 0:
        return va_corr

    # NOTE: it is assumed that v2, v3 angles and va scale are small enough
    # so that for expected scale factors the issue of angle wrapping
    # (180 degrees) can be neglected.
    v2_shift = (1 - va_scale) * v2_ref
    v3_shift = (1 - va_scale) * v3_ref

    va_corr |= Shift(v2_shift, name="dva_v2_shift") & Shift(
        v3_shift, name="dva_v3_shift"
    )
    va_corr.name = "DVA_Correction"
    return va_corr
