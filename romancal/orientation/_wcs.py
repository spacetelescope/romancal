"""World Coordinate System (WCS) utilities"""

import logging

import numpy as np
from astropy.time import Time
from stcal.alignment.util import compute_s_region_keyword

from . import _lib as olib
from . import _pointing as plib
from . import _siaf as siaf_lib
from . import _transforms as tlib

__all__ = []

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LOGLEVELS = [logging.INFO, logging.DEBUG, olib.DEBUG_FULL]


def calc_wcs(t_pars: tlib.TransformParameters):
    """
    Calculate WCS.

    Given observatory orientation and target aperture,
    calculate V1 and Reference Pixel sky coordinates.

    Parameters
    ----------
    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    wcsinfo, vinfo, transforms : WCSRef, WCSRef, Transforms
        A 3-tuple is returned with the WCS pointing for
        the aperture and the V1 axis, and the transformation matrices.
    """
    # Calculate transforms
    transforms = tlib.calc_transforms(t_pars)

    # Calculate the V1 WCS information
    vinfo = tlib.calc_wcs_from_matrix(transforms.m_eci2v)

    # Calculate the Aperture WCS
    wcsinfo = wcsinfo_from_siaf(t_pars.aperture, vinfo)

    # That's all folks
    return wcsinfo, vinfo, transforms


def wcsinfo_from_siaf(aperture, vinfo):
    """Calculate aperture reference point WCS from V-frame WCS and SIAF

    Parameters
    ----------
    aperture : str
        The aperture in use

    vinfo : WCSRef
        The V-frame WCS

    Returns
    -------
    wcsinfo : WCSRef
        The WCS for the aperture's reference point, as defined by its SIAF.
    """
    from pysiaf.utils.rotations import sky_posangle

    wfi = siaf_lib.SIAF[aperture.upper()]

    # For transformations between the telescope frame and all other frames,
    # an attitude matrix is created using the V-frame WCS information.
    attitude = olib.attitude_from_v1(vinfo)
    wfi.set_attitude_matrix(attitude)
    skycoord = wfi.reference_point(to_frame="sky")
    pa_v3 = sky_posangle(attitude, *skycoord)

    # Compute S_REGION keyword
    corners = wfi.corners("sky")
    footprint = np.array([corners[0], corners[1]]).T
    s_region_keyword = compute_s_region_keyword(footprint)

    wcsinfo = olib.WCSRef(
        ra=skycoord[0], dec=skycoord[1], pa=pa_v3, s_region=s_region_keyword
    )
    return wcsinfo


def update_wcs_from_telem(model, t_pars: tlib.TransformParameters):
    """
    Update WCS pointing information.

    Given a `roman.datamodels.DataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present.

    Parameters
    ----------
    model : `~roman.datamodels.DataModel`
        The model to update. The update is done in-place.

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms or None
        If available, the transformation matrices.
    """
    logger.info("Updating wcs from telemetry.")

    # Initialization. If provided, provide a default Pointing.
    transforms = None  # Assume no transforms are calculated.
    quality = None  # Unknown pointing quality.

    # Setup SIAF information.
    if t_pars.siaf_path is not None:
        logger.info("Using pysiaf xml folder %s", t_pars.siaf_path)
    else:
        logger.info("Using build-in pysiaf xml")
    siaf_lib.open_siaf(basepath=t_pars.siaf_path)

    # Get the pointing information
    try:
        t_pars.update_from_engdb()
    except ValueError as exception:
        logger.error("Cannot retrieve valid engineering orientation data")
        if t_pars.default_quaternion is None or not t_pars.allow_default:
            logger.error("Use of default orientation has been disabled. Aborting.")
            raise
        else:
            logger.warning("Exception is %s", exception)
            obstime = Time(
                (t_pars.obsstart.mjd + t_pars.obsend.mjd) / 2.0, format="mjd"
            )
            logger.warning(
                "Using provided default quaternion: %s", t_pars.default_quaternion
            )
            logger.warning("    at time %s", obstime.iso)
            logger.info("Setting pointing quality to PLANNED")
            t_pars.pointing = plib.Pointing(
                q=t_pars.default_quaternion, obstime=obstime
            )
            quality = "PLANNED"
    else:
        logger.info("Successful read of engineering quaternions:")
        logger.info("\tPointing: %s", t_pars.pointing)
        quality = "CALCULATED"

    # Attempt to calculate WCS information
    try:
        wcsinfo, vinfo, transforms = calc_wcs(t_pars)
        logger.info("Setting pointing quality to %s", quality)
    except olib.EXPECTED_ERRORS:
        logger.error("WCS calculation has failed")
        raise

    # Update model meta.
    logger.info("Aperture WCS info: %s", wcsinfo)
    logger.info("V1 WCS info: %s", vinfo)
    olib.update_meta(model, t_pars, wcsinfo, vinfo, quality)

    return transforms


def calc_wcs_over_time(obsstart, obsend, t_pars: tlib.TransformParameters):
    """
    Calculate V1 and aperture WCS over a time period.

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    obstimes, wcsinfos, vinfos : [astropy.time.Time[,...]], [WCSRef[,...]], [WCSRef[,...]]
        A 3-tuple is returned with the WCS pointings for
        the aperture and the V1 axis.
    """
    # Setup structures
    obstimes = []
    wcsinfos = []
    vinfos = []

    # Calculate WCS
    try:
        pointings = plib.get_pointing(
            obsstart,
            obsend,
            service_kwargs=t_pars.service_kwargs,
            tolerance=t_pars.tolerance,
            reduce_func=t_pars.reduce_func,
        )
    except ValueError:
        logger.warning(
            "Cannot get valid engineering mnemonics from engineering database"
        )
        raise
    if not isinstance(pointings, list):
        pointings = [pointings]
    for pointing in pointings:
        t_pars.pointing = pointing
        wcsinfo, vinfo, transforms = calc_wcs(t_pars)
        obstimes.append(pointing.obstime)
        wcsinfos.append(wcsinfo)
        vinfos.append(vinfo)

    return obstimes, wcsinfos, vinfos
