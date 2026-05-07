"""
WCS construction utilities for Roman WFI images.
"""

import logging

import gwcs.coordinate_frames as cf
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import bind_bounding_box
from astropy.modeling.models import Identity, RotationSequence3D, Scale, Shift
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from gwcs.wcs import WCS, Step
from roman_datamodels import datamodels as rdm
from stcal.alignment.util import compute_s_region_keyword, wcs_bbox_from_shape

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def load_wcs(input_model, reference_files=None):
    """Create a gWCS object and store it in ``Model.meta``.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.WfiImage`
        The exposure.
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.

    Returns
    -------
    output_model : `~roman_datamodels.ImageModel`
        The input image file with attached gWCS object.
        The input_model is modified in place.
    """
    output_model = input_model

    if reference_files is not None:
        for ref_type, ref_file in reference_files.items():
            reference_files[ref_type] = (
                ref_file if ref_file not in ["N/A", ""] else None
            )
    else:
        reference_files = {}

    # Frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    v2v3vacorr = cf.Frame2D(
        name="v2v3vacorr",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # Transforms between frames
    distortion = _wfi_distortion(output_model, reference_files)
    tel2sky = v23tosky(output_model)

    # Compute differential velocity aberration (DVA) correction:
    va_corr = _dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )

    pipeline = [
        Step(detector, distortion),
        Step(v2v3, va_corr),
        Step(v2v3vacorr, tel2sky),
        Step(world, None),
    ]
    wcs = WCS(pipeline)
    if wcs.bounding_box is None:
        wcs.bounding_box = wcs_bbox_from_shape(output_model.data.shape)

    output_model.meta["wcs"] = wcs

    # update S_REGION
    add_s_region(output_model)

    output_model.meta.cal_step["assign_wcs"] = "COMPLETE"

    return output_model


def _wfi_distortion(model, reference_files):
    """
    Create the "detector" to "v2v3" transform for WFI

    Parameters
    ----------
    model : `~roman_datamodels.datamodels.WfiImage`
        The data model for processing
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.

    Returns
    -------
    The transform model
    """

    dist = rdm.DistortionRefModel(reference_files["distortion"])
    transform = dist.coordinate_distortion_transform

    try:
        bbox = transform.bounding_box.bounding_box(order="F")
    except NotImplementedError:
        # Check if the transform in the reference file has a ``bounding_box``.
        # If not set a ``bounding_box`` equal to the size of the image after
        # assembling all distortion corrections.
        bbox = None
    dist.close()

    bind_bounding_box(
        transform,
        wcs_bbox_from_shape(model.data.shape) if bbox is None else bbox,
        order="F",
    )

    return transform


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


def _dva_corr_model(va_scale, v2_ref, v3_ref):
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


def _create_footprint(wcs, shape=None, center=False):
    """Calculate sky footprint

    Parameters
    ----------
    wcs : `gwcs.WCS`
        The WCS information to get the footprint from

    shape : n-tuple or None
       Shape to use if wcs has no defined shape.

    center : bool
        If True use the center of the pixel, otherwise use the corner.

    Returns
    -------
    footprint : `numpy.ndarray`
        The footprint.
    """
    bbox = wcs.bounding_box

    if bbox is None:
        bbox = wcs_bbox_from_shape(shape)

    # footprint is an array of shape (2, 4) - i.e. 4 values for RA and 4 values for
    # Dec - as we are interested only in the footprint on the sky
    footprint = wcs.footprint(bbox, center=center, axis_type="spatial").T
    # take only imaging footprint
    footprint = footprint[:2, :]

    # Make sure RA values are all positive
    negative_ind = footprint[0] < 0
    if negative_ind.any():
        footprint[0][negative_ind] = 360 + footprint[0][negative_ind]

    footprint = footprint.T
    return footprint


def add_s_region(model):
    """
    Calculate the detector's footprint using ``WCS.footprint`` and save it in the
    ``S_REGION`` keyword

    Parameters
    ----------
    model : `~roman_datamodels.datamodels.ImageModel`
        The data model for processing

    Returns
    -------
    A formatted string representing the detector's footprint
    """
    _update_s_region_keyword(
        model, _create_footprint(model.meta.wcs, shape=model.shape, center=False)
    )


def _update_s_region_keyword(model, footprint):
    s_region = compute_s_region_keyword(footprint)
    log.info(f"S_REGION VALUES: {s_region}")
    if "nan" in s_region:
        # do not update s_region if there are NaNs.
        log.info("There are NaNs in s_region, S_REGION not updated.")
    else:
        model.meta.wcsinfo.s_region = s_region
        log.info(f"Update S_REGION to {model.meta.wcsinfo.s_region}")
