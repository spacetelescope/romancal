import functools
import logging
from typing import List, Tuple, Union

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.modeling import models as astmodels
from astropy.utils.misc import isiterable
from gwcs import WCS
from gwcs.wcstools import wcs_from_fiducial
from roman_datamodels.datamodels import DataModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_MAX_SIP_DEGREE = 6


def wcs_bbox_from_shape(shape):
    """Create a bounding box from the shape of the data.
    This is appropriate to attach to a wcs object.

    Parameters
    ----------
    shape : tuple
        The shape attribute from a `numpy.ndarray` array.

    Returns
    -------
    bbox : tuple
        Bounding box in x, y order.
    """
    bbox = ((-0.5, shape[-1] - 0.5), (-0.5, shape[-2] - 0.5))
    return bbox


def wcs_from_footprints(
    dmodels,
    refmodel=None,
    transform=None,
    bounding_box=None,
    pscale_ratio=None,
    pscale=None,
    rotation=None,
    shape=None,
    ref_pixel: Tuple[float, float] = None,
    ref_coord: Tuple[float, float] = None,
):
    """
    Create a WCS from a list of input data models.

    A fiducial point in the output coordinate frame is created from  the
    footprints of all WCS objects. For a spatial frame this is the center
    of the union of the footprints. For a spectral frame the fiducial is in
    the beginning of the footprint range.
    If ``refmodel`` is None, the first WCS object in the list is considered
    a reference. The output coordinate frame and projection (for celestial frames)
    is taken from ``refmodel``.
    If ``transform`` is not supplied, a compound transform is created using
    CDELTs and PC.
    If ``bounding_box`` is not supplied, the bounding_box of the new WCS is computed
    from bounding_box of all input WCSs.

    Parameters
    ----------
    dmodels : list of `~jwst.datamodels.DataModel`
        A list of data models.
    refmodel : `~jwst.datamodels.DataModel`, optional
        This model's WCS is used as a reference.
        WCS. The output coordinate frame, the projection and a
        scaling and rotation transform is created from it. If not supplied
        the first model in the list is used as ``refmodel``.
    transform : `~astropy.modeling.core.Model`, optional
        A transform, passed to :meth:`~gwcs.wcstools.wcs_from_fiducial`
        If not supplied Scaling | Rotation is computed from ``refmodel``.
    bounding_box : tuple, optional
        Bounding_box of the new WCS.
        If not supplied it is computed from the bounding_box of all inputs.
    pscale_ratio : float, None, optional
        Ratio of input to output pixel scale. Ignored when either
        ``transform`` or ``pscale`` are provided.
    pscale : float, None, optional
        Absolute pixel scale in degrees. When provided, overrides
        ``pscale_ratio``. Ignored when ``transform`` is provided.
    rotation : float, None, optional
        Position angle (in degrees) of output image's Y-axis relative to North.
        A value of 0.0 would orient the final output image to be North up.
        The default of `None` specifies that the images will not be rotated,
        but will instead be resampled in the default orientation for the camera
        with the x and y axes of the resampled image corresponding
        approximately to the detector axes. Ignored when ``transform`` is
        provided.
    shape : tuple of int, None, optional
        Shape of the image (data array) using ``numpy.ndarray`` convention
        (``ny`` first and ``nx`` second). This value will be assigned to
        ``pixel_shape`` and ``array_shape`` properties of the returned
        WCS object.
    ref_pixel : tuple of float, None, optional
        Position of the reference pixel in the image array.  If ``ref_pixel`` is not
        specified, it will be set to the center of the bounding box of the
        returned WCS object.
    ref_coord : tuple of float, None, optional
        Right ascension and declination of the reference pixel. Automatically
        computed if not provided.

    """
    bb = bounding_box
    wcslist = [im.meta.wcs for im in dmodels]

    if not isiterable(wcslist):
        raise ValueError("Expected 'wcslist' to be an iterable of WCS objects.")

    if not all([isinstance(w, WCS) for w in wcslist]):
        raise TypeError("All items in wcslist are to be instances of gwcs.WCS.")

    if refmodel is None:
        refmodel = dmodels[0]
    else:
        if not isinstance(refmodel, DataModel):
            raise TypeError("Expected refmodel to be an instance of DataModel.")

    fiducial = compute_fiducial(wcslist, bb)
    if ref_coord is not None:
        # overwrite spatial axes with user-provided ref_coord:
        i = 0
        for k, axt in enumerate(wcslist[0].output_frame.axes_type):
            if axt == "SPATIAL":
                fiducial[k] = ref_coord[i]
                i += 1

    ref_fiducial = np.array(
        [refmodel.meta.wcsinfo.ra_ref, refmodel.meta.wcsinfo.dec_ref]
    )

    prj = astmodels.Pix2Sky_TAN()

    if transform is None:
        transform = []
        sky_axes = refmodel.meta.wcs._get_axes_indices().tolist()

        # Need to put the rotation matrix (List[float, float, float, float])
        # returned from calc_rotation_matrix into the correct shape for
        # constructing the transformation
        v3yangle = np.deg2rad(refmodel.meta.wcsinfo.v3yangle)
        vparity = refmodel.meta.wcsinfo.vparity
        if rotation is None:
            roll_ref = np.deg2rad(refmodel.meta.wcsinfo.roll_ref)
        else:
            roll_ref = np.deg2rad(rotation) + (vparity * v3yangle)

        pc = np.reshape(
            calc_rotation_matrix(roll_ref, v3yangle, vparity=vparity), (2, 2)
        )

        rotation = astmodels.AffineTransformation2D(pc, name="pc_rotation_matrix")
        transform.append(rotation)

        if sky_axes:
            if not pscale:
                pscale = compute_scale(
                    refmodel.meta.wcs, ref_fiducial, pscale_ratio=pscale_ratio
                )
            transform.append(
                astmodels.Scale(pscale, name="cdelt1")
                & astmodels.Scale(pscale, name="cdelt2")
            )

        if transform:
            transform = functools.reduce(lambda x, y: x | y, transform)

    out_frame = refmodel.meta.wcs.output_frame
    input_frame = refmodel.meta.wcs.input_frame
    wnew = wcs_from_fiducial(
        fiducial,
        coordinate_frame=out_frame,
        projection=prj,
        transform=transform,
        input_frame=input_frame,
    )

    footprints = [w.footprint().T for w in wcslist]
    domain_bounds = np.hstack([wnew.backward_transform(*f) for f in footprints])
    axis_min_values = np.min(domain_bounds, axis=1)
    domain_bounds = (domain_bounds.T - axis_min_values).T

    output_bounding_box = []
    for axis in out_frame.axes_order:
        axis_min, axis_max = (
            domain_bounds[axis].min(),
            domain_bounds[axis].max(),
        )
        output_bounding_box.append((axis_min, axis_max))

    output_bounding_box = tuple(output_bounding_box)
    if ref_pixel is None:
        offset1, offset2 = wnew.backward_transform(*fiducial)
        offset1 -= axis_min_values[0]
        offset2 -= axis_min_values[1]
    else:
        offset1, offset2 = ref_pixel
    offsets = astmodels.Shift(-offset1, name="ref_pixel1") & astmodels.Shift(
        -offset2, name="ref_pixel2"
    )

    wnew.insert_transform("detector", offsets, after=True)
    wnew.bounding_box = output_bounding_box

    if shape is None:
        shape = [int(axs[1] - axs[0] + 0.5) for axs in output_bounding_box[::-1]]

    wnew.pixel_shape = shape[::-1]
    wnew.array_shape = shape

    return wnew


def compute_scale(
    wcs: WCS,
    fiducial: Union[tuple, np.ndarray],
    disp_axis: int = None,
    pscale_ratio: float = None,
) -> float:
    """Compute scaling transform.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        Reference WCS object from which to compute a scaling factor.

    fiducial : tuple
        Input fiducial of (RA, DEC) or (RA, DEC, Wavelength) used in calculating
        reference points.

    disp_axis : int
        Dispersion axis integer. Assumes the same convention as
        `wcsinfo.dispersion_direction`

    pscale_ratio : int
        Ratio of input to output pixel scale

    Returns
    -------
    scale : float
        Scaling factor for x and y or cross-dispersion direction.

    """
    spectral = "SPECTRAL" in wcs.output_frame.axes_type

    if spectral and disp_axis is None:
        raise ValueError("If input WCS is spectral, a disp_axis must be given")

    crpix = np.array(wcs.invert(*fiducial))

    delta = np.zeros_like(crpix)
    spatial_idx = np.where(np.array(wcs.output_frame.axes_type) == "SPATIAL")[0]
    delta[spatial_idx[0]] = 1

    crpix_with_offsets = np.vstack((crpix, crpix + delta, crpix + np.roll(delta, 1))).T
    crval_with_offsets = wcs(*crpix_with_offsets, with_bounding_box=False)

    coords = SkyCoord(
        ra=crval_with_offsets[spatial_idx[0]],
        dec=crval_with_offsets[spatial_idx[1]],
        unit="deg",
    )
    xscale = np.abs(coords[0].separation(coords[1]).value)
    yscale = np.abs(coords[0].separation(coords[2]).value)

    if pscale_ratio is not None:
        xscale *= pscale_ratio
        yscale *= pscale_ratio

    if spectral:
        # Assuming scale doesn't change with wavelength
        # Assuming disp_axis is consistent with
        # DataModel.meta.wcsinfo.dispersion.direction
        return yscale if disp_axis == 1 else xscale

    return np.sqrt(xscale * yscale)


def calc_rotation_matrix(
    roll_ref: float, v3i_yang: float, vparity: int = 1
) -> List[float]:
    """Calculate the rotation matrix.

    Parameters
    ----------
    roll_ref : float
        Telescope roll angle of V3 North over East at the ref. point in radians

    v3i_yang : float
        The angle between ideal Y-axis and V3 in radians.

    vparity : int
        The x-axis parity, usually taken from the JWST SIAF parameter VIdlParity.
        Value should be "1" or "-1".

    Returns
    -------
    matrix: [pc1_1, pc1_2, pc2_1, pc2_2]
        The rotation matrix

    Notes
    -----
    The rotation is

       ----------------
       | pc1_1  pc2_1 |
       | pc1_2  pc2_2 |
       ----------------

    """
    if vparity not in (1, -1):
        raise ValueError(f"vparity should be 1 or -1. Input was: {vparity}")

    rel_angle = roll_ref - (vparity * v3i_yang)

    pc1_1 = vparity * np.cos(rel_angle)
    pc1_2 = np.sin(rel_angle)
    pc2_1 = vparity * -np.sin(rel_angle)
    pc2_2 = np.cos(rel_angle)

    return [pc1_1, pc1_2, pc2_1, pc2_2]


def compute_fiducial(wcslist, bounding_box=None):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.
    """

    axes_types = wcslist[0].output_frame.axes_type
    spatial_axes = np.array(axes_types) == "SPATIAL"
    spectral_axes = np.array(axes_types) == "SPECTRAL"
    footprints = np.hstack([w.footprint(bounding_box=bounding_box).T for w in wcslist])
    spatial_footprint = footprints[spatial_axes]
    spectral_footprint = footprints[spectral_axes]

    fiducial = np.empty(len(axes_types))
    if spatial_footprint.any():
        lon, lat = spatial_footprint
        lon, lat = np.deg2rad(lon), np.deg2rad(lat)
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)

        x_mid = (np.max(x) + np.min(x)) / 2.0
        y_mid = (np.max(y) + np.min(y)) / 2.0
        z_mid = (np.max(z) + np.min(z)) / 2.0
        lon_fiducial = np.rad2deg(np.arctan2(y_mid, x_mid)) % 360.0
        lat_fiducial = np.rad2deg(np.arctan2(z_mid, np.sqrt(x_mid**2 + y_mid**2)))
        fiducial[spatial_axes] = lon_fiducial, lat_fiducial
    if spectral_footprint.any():
        fiducial[spectral_axes] = spectral_footprint.min()
    return fiducial


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

    bbox = model.meta.wcs.bounding_box

    if bbox is None:
        bbox = wcs_bbox_from_shape(model.data.shape)

    # footprint is an array of shape (2, 4) - i.e. 4 values for RA and 4 values for
    # Dec - as we are interested only in the footprint on the sky
    footprint = model.meta.wcs.footprint(bbox, center=True, axis_type="spatial").T
    # take only imaging footprint
    footprint = footprint[:2, :]

    # Make sure RA values are all positive
    negative_ind = footprint[0] < 0
    if negative_ind.any():
        footprint[0][negative_ind] = 360 + footprint[0][negative_ind]

    footprint = footprint.T
    update_s_region_keyword(model, footprint)


def update_s_region_keyword(model, footprint):
    s_region = "POLYGON ICRS " + " ".join([str(x) for x in footprint.ravel()]) + " "
    log.info(f"S_REGION VALUES: {s_region}")
    if "nan" in s_region:
        # do not update s_region if there are NaNs.
        log.info("There are NaNs in s_region, S_REGION not updated.")
    else:
        model.meta.wcsinfo.s_region = s_region
        log.info(f"Update S_REGION to {model.meta.wcsinfo.s_region}")
