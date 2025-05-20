import logging

import numpy as np
from astropy.coordinates import SkyCoord
from gwcs import WCS
from stcal.alignment.util import compute_s_region_keyword

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


def compute_scale(
    wcs: WCS,
    fiducial: tuple | np.ndarray,
    disp_axis: int | None = None,
    pscale_ratio: float | None = None,
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
        Ratio of output pixel scale to input pixel scale.

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
) -> list[float]:
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


def create_footprint(wcs, shape=None, center=False):
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
    update_s_region_keyword(
        model, create_footprint(model.meta.wcs, shape=model.shape, center=False)
    )


def update_s_region_keyword(model, footprint):
    s_region = compute_s_region_keyword(footprint)
    log.info(f"S_REGION VALUES: {s_region}")
    if "nan" in s_region:
        # do not update s_region if there are NaNs.
        log.info("There are NaNs in s_region, S_REGION not updated.")
    else:
        model.meta.wcsinfo.s_region = s_region
        log.info(f"Update S_REGION to {model.meta.wcsinfo.s_region}")


def list_1d_to_2d(l, n):
    """Convert 1-dimensional list to 2-dimensional

    Parameters
    ----------
    l : list
        The list to convert.

    n : int
       The length of the x dimension, or the length of the inner lists.

    Returns
    -------
    l2d : list of lists
        The 2D form
    """
    return [l[i : i + n] for i in range(0, len(l), n)]
