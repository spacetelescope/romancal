"""World Coordinate System (WCS) utilities"""

import logging

import asdf
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.time import Time
from astropy import units as u
import numpy as np
import roman_datamodels as rdm
from stcal.alignment.util import compute_s_region_keyword
from scipy.spatial.transform import Rotation as R

from . import _lib as olib
from . import _pointing as plib
from . import _siaf as siaf_lib
from . import _transforms as tlib

__all__ = []

# M_XtoZ matrix for the forward construction: M_OBStoFCS = M_XtoZ @ Rx(ya) @ Ry(bz) @ Rz(-by)
M_X2Z = np.array(
    [
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
    ]
)

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LOGLEVELS = [logging.INFO, logging.DEBUG, olib.DEBUG_FULL]


def apply_vector_velocity_aberration(u_vec, v_xyz_kms):
    """Applies velocity aberration to a 3D ECI unit vector."""
    c_kms = 299792.458
    beta = np.array(v_xyz_kms, dtype=float) / c_kms

    beta_dot_u = np.dot(beta, u_vec)
    u_app = u_vec + beta - (beta_dot_u * u_vec)
    u_app /= np.linalg.norm(u_app)

    scale_factor = 1.0 + beta_dot_u
    return u_app, float(scale_factor)

def build_v2v3_skyoffset_frame(v1_vec, v3_vec):
    """
    Constructs an astropy SkyOffsetFrame from V1 and V3 ECI unit vectors,
    returning (v2v3_frame, ra_v1_deg, dec_v1_deg, pa_v3_deg).
    """
    ra_v1 = float(np.rad2deg(np.arctan2(v1_vec[1], v1_vec[0])) % 360.0)
    dec_v1 = float(np.rad2deg(np.arcsin(np.clip(v1_vec[2], -1.0, 1.0))))

    # Infinitesimal step along +V3 to measure PA_V3 East of North via astropy
    v1_coord = SkyCoord(ra=ra_v1 * u.deg, dec=dec_v1 * u.deg, frame="icrs")
    v3_step = v1_vec + 1e-6 * v3_vec
    v3_step /= np.linalg.norm(v3_step)
    v3_step_coord = SkyCoord(
        ra=np.rad2deg(np.arctan2(v3_step[1], v3_step[0])) * u.deg,
        dec=np.rad2deg(np.arcsin(np.clip(v3_step[2], -1.0, 1.0))) * u.deg,
        frame="icrs",
    )
    pa_v3 = float(v1_coord.position_angle(v3_step_coord).deg % 360.0)

    # SkyOffsetFrame positive rotation is CCW (West of North).
    # Passing +pa_v3 aligns +lat with +V3 (East of North) and +lon with +V2.
    v2v3_frame = SkyOffsetFrame(origin=v1_coord, rotation=pa_v3 * u.deg)
    return v2v3_frame, ra_v1, dec_v1, pa_v3


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
    meta : dict
    Set of meta data that defines the wcs. The primary keys are:
    - pointing: The pointing meta block
    - velocity_aberration: The aberration characteristics
    - wcsinfo: The basic WCS information

    """
    # 0. LOAD REFERENCE FILES
    with asdf.open(t_pars.bam_ref) as af_bam:
        q_bam = [float(x) for x in af_bam["roman"]["bam_quaternion"]]

    sca_info = extract_sca_wcs_params(t_pars.distortion_ref)
    aperture_name = t_pars.aperture
    if aperture_name is None:
        aperture_name = f"{sca_info['detector']}_FULL"

    # 1. APPARENT ATTITUDE FRAME (From Telemetry qbj)
    # The physical orientation of the telescope capturing aberrated photons
    att_matrix_app = R.from_quat(t_pars.pointing.q).as_matrix()
    v1_app = att_matrix_app[:, 0]
    v3_app = att_matrix_app[:, 2]

    app_v2v3_frame, _, _, _ = build_v2v3_skyoffset_frame(v1_app, v3_app)

    # 2. TRUE ATTITUDE FRAME (Inverse Velocity Aberration)
    # Apply -v_xyz_kms to back out the aberration and recover the ICRS catalog frame
    v_xyz_kms_inv = -np.array(t_pars.velocity, dtype=float)
    v1_true, _ = apply_vector_velocity_aberration(v1_app, v_xyz_kms_inv)
    v3_true_approx, _ = apply_vector_velocity_aberration(v3_app, v_xyz_kms_inv)

    # Re-orthogonalize V3 against V1 to enforce a rigid 90-degree body frame
    v3_true = v3_true_approx - np.dot(v3_true_approx, v1_true) * v1_true
    v3_true /= np.linalg.norm(v3_true)

    true_v2v3_frame, true_ra_v1, true_dec_v1, true_pa_v3 = build_v2v3_skyoffset_frame(
        v1_true, v3_true
    )

    # 3. GLOBAL TARGET POINTING VIA BAM QUATERNION
    tgt_ya_deg, tgt_v2_ref, tgt_v3_ref = get_angles_from_quaternion(q_bam)

    # Project True Sky Coordinates for APT alignment
    tgt_true_ra, tgt_true_dec, true_roll_target, tgt_true_coord = (
        project_v2v3_with_roll(tgt_v2_ref, tgt_v3_ref, true_v2v3_frame)
    )
    pa_aperture_true = (true_roll_target + tgt_ya_deg) % 360.0

    # Velocity Aberration Scale Factor using True Target SkyCoord (forward velocity)
    u_tgt_true = tgt_true_coord.cartesian.xyz.value
    _, global_scale_factor = apply_vector_velocity_aberration(u_tgt_true, t_pars.velocity)

    # Apparent Reference Coordinates at Target
    tgt_app_coord = SkyCoord(
        lon=tgt_v2_ref * u.arcsec, lat=tgt_v3_ref * u.arcsec, frame=app_v2v3_frame
    ).transform_to("icrs")

    # 4. SCA-SPECIFIC WCS PROJECTIONS VIA DISTORTION REFERENCE FILE (True Frame)
    true_ra_ref, true_dec_ref, true_roll_ref, _ = project_v2v3_with_roll(
        sca_info["v2_ref"], sca_info["v3_ref"], true_v2v3_frame
    )

    corners_coord = SkyCoord(
        lon=sca_info["corners_v2"] * u.arcsec,
        lat=sca_info["corners_v3"] * u.arcsec,
        frame=true_v2v3_frame,
    ).transform_to("icrs")

    polygon_str = " ".join(
        [
            f"{ra:.9f} {dec:.9f}"
            for ra, dec in zip(corners_coord.ra.deg, corners_coord.dec.deg, strict=True)
        ]
    )

    # 5. ASDF METADATA ASSEMBLY
    wcs_meta = {
        "pointing": {
            "pointing_engineering_source": "CALCULATED",
            "quaternion": [float(x) for x in t_pars.pointing.q],
            "bam_quaternion": q_bam,
            "ra_v1": true_ra_v1,
            "dec_v1": true_dec_v1,
            "pa_v3": true_pa_v3,
            "target_aperture": t_pars.target_aperture,
            "target_ra": tgt_true_ra,
            "target_dec": tgt_true_dec,
            "pa_aperture": float(pa_aperture_true),
        },
        "velocity_aberration": {
            "ra_reference": float(tgt_app_coord.ra.deg),
            "dec_reference": float(tgt_app_coord.dec.deg),
            "scale_factor": global_scale_factor,
        },
        "wcsinfo": {
            "aperture_name": aperture_name,
            "v2_ref": sca_info["v2_ref"],
            "v3_ref": sca_info["v3_ref"],
            "v3yangle": sca_info["v3yangle"],
            "vparity": sca_info["vparity"],
            "ra_ref": true_ra_ref,
            "dec_ref": true_dec_ref,
            "roll_ref": true_roll_ref,
            "s_region": f"POLYGON ICRS {polygon_str}",
        },
    }
    return wcs_meta

def get_angles_from_quaternion(q_inst2b):
    """
    Recover Roman focal-plane geometry angles from an Instrument -> Body quaternion.
    Returns V3IdlYAngle (deg), V2Ref (arcsec), and V3Ref (arcsec).
    """
    q_inst2b = np.asarray(q_inst2b, dtype=float)
    q_inst2b /= np.linalg.norm(q_inst2b)

    m_b2inst = R.from_quat(q_inst2b).inv().as_matrix()
    a = M_X2Z.T @ m_b2inst

    bz_rad = np.arcsin(np.clip(a[0, 2], -1.0, 1.0))
    by_rad = np.arctan2(-a[0, 1], a[0, 0])
    ya_rad = np.arctan2(-a[1, 2], a[2, 2])

    ya_deg = np.rad2deg(ya_rad)
    by_arcsec = np.rad2deg(by_rad) * 3600.0
    bz_arcsec = np.rad2deg(bz_rad) * 3600.0

    if ya_deg > 90.0:
        ya_deg -= 180.0

    ya_deg = (ya_deg + 180.0) % 360.0 - 180.0
    by_arcsec *= -1.0

    return float(ya_deg), float(by_arcsec), float(bz_arcsec)


def project_v2v3_with_roll(v2_arcsec, v3_arcsec, v2v3_frame):
    """
    Transforms (V2, V3) in arcsec to ICRS (RA, Dec) and computes the local V3
    position angle (roll) using astropy SkyCoord.
    """
    c_ref = SkyCoord(
        lon=v2_arcsec * u.arcsec, lat=v3_arcsec * u.arcsec, frame=v2v3_frame
    )
    c_icrs = c_ref.transform_to("icrs")

    # Infinitesimal step along local +V3 meridian (+lat in SkyOffsetFrame)
    c_north_v3 = SkyCoord(
        lon=v2_arcsec * u.arcsec,
        lat=(v3_arcsec + 1e-3) * u.arcsec,
        frame=v2v3_frame,
    ).transform_to("icrs")

    local_roll = float(c_icrs.position_angle(c_north_v3).deg % 360.0)
    return float(c_icrs.ra.deg), float(c_icrs.dec.deg), local_roll, c_icrs


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
    wcs_meta : dict
        The various wcs-related meta keywords calculated for.
    """
    logger.info("Updating wcs from telemetry.")

    # Unknown pointing quality.
    quality = None

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
        wcs_meta = calc_wcs(t_pars)
        logger.info("Setting pointing quality to %s", quality)
    except olib.EXPECTED_ERRORS:
        logger.error("WCS calculation has failed")
        raise

    # Update model meta.
    logger.info("Updated WCS meta information: %s", wcs_meta)
    olib.update_meta(model, t_pars, wcs_meta, quality)

    return wcs_meta


# ###
# old code
# ###
def wcsinfo_from_siaf(aperture, vinfo, siaf_path=None):
    """Calculate aperture reference point WCS from V-frame WCS and SIAF

    Parameters
    ----------
    aperture : str
        The aperture in use

    vinfo : WCSRef
        The V-frame WCS

    siaf_path : str
        The folder where the SIAF xml information resides.

    Returns
    -------
    wcsinfo : WCSRef
        The WCS for the aperture's reference point, as defined by its SIAF.
    """
    from pysiaf.utils.rotations import sky_posangle

    siaf = siaf_lib.open_siaf(basepath=siaf_path)
    wfi = siaf[aperture.upper()]

    # For transformations between the telescope frame and all other frames,
    # an attitude matrix is created using the V-frame WCS information.
    attitude = olib.attitude_from_v1(vinfo, siaf_path=siaf_path)
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


def extract_sca_wcs_params(distortion_path):
    """
    Extracts V2Ref, V3Ref, V3IdlYAngle, VIdlParity, detector name, and the
    4 corner (V2, V3) coordinates using roman_datamodels.
    """
    corners_x_0idx = np.array([3.5, 4091.5, 4091.5, 3.5])
    corners_y_0idx = np.array([3.5, 3.5, 4091.5, 4091.5])

    with rdm.open(distortion_path) as dm:
        detector = str(dm.meta.instrument.detector)
        dist_mod = dm.coordinate_distortion_transform

        # Direct parameter access on the CompoundModel
        v2_ref = float(dist_mod.offset_10.value)
        v3_ref = float(dist_mod.offset_11.value)

        # Ideal -> V2/V3 rotation & parity parameters (submodels 8 and 9)
        c10_v2 = float(dist_mod.c1_0_8.value)  # p * cos(ya)
        c01_v2 = float(dist_mod.c0_1_8.value)  # sin(ya)
        c10_v3 = float(dist_mod.c1_0_9.value)  # -p * sin(ya)
        c01_v3 = float(dist_mod.c0_1_9.value)  # cos(ya)

        v3_idl_yangle = float(np.rad2deg(np.arctan2(c01_v2, c01_v3)))
        det_rot = (c10_v2 * c01_v3) - (c01_v2 * c10_v3)
        v_idl_parity = int(np.round(det_rot))

        corners_v2, corners_v3 = dist_mod(corners_x_0idx, corners_y_0idx)

        return {
            "detector": detector,
            "v2_ref": v2_ref,
            "v3_ref": v3_ref,
            "v3yangle": v3_idl_yangle,
            "vparity": v_idl_parity,
            "corners_v2": np.array(corners_v2, dtype=float),
            "corners_v3": np.array(corners_v3, dtype=float),
        }
