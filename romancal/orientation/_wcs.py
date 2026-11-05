"""World Coordinate System (WCS) utilities"""

import logging

import asdf
import numpy as np
import roman_datamodels as rdm
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.time import Time
from scipy.spatial.transform import Rotation as R
from stcal.velocity_aberration import compute_va_effects_vector

from . import _lib as olib
from . import _pointing as plib
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
    wcs_meta : dict
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

    # Check that observatory velocity is "reasonable" for a given definition of reasonable.
    velocity = np.array(t_pars.velocity)
    if (
        None in velocity
        or np.all(velocity == 0.0)
        or (np.abs(velocity) > olib.MAX_OBSERVATORY_SPEED).any()
    ):
        logger.warning(
            "Velocity: %s is either unspecified or contains unreasonable values. Assuming zero velocity.",
            velocity,
        )
        velocity = np.zeros(velocity.shape)

    # 1. APPARENT ATTITUDE FRAME (From Telemetry qbj)
    att_matrix_app = R.from_quat(t_pars.pointing.q).as_matrix()
    v1_app = att_matrix_app[:, 0]
    v3_app = att_matrix_app[:, 2]

    app_v2v3_frame, _, _, _ = build_v2v3_skyoffset_frame(v1_app, v3_app)

    # 2. TRUE (UN-ABERRATED) V1 BORESIGHT & FRAME
    def true_point(v2, v3):
        return v2v3_to_true_icrs(v2, v3, app_v2v3_frame, velocity)

    true_ra_v1, true_dec_v1, true_pa_v3, _ = true_point_and_roll(
        0.0, 0.0, app_v2v3_frame, velocity
    )

    # 3. GLOBAL TARGET POINTING (DYNAMIC APERTURE RESOLUTION)
    if t_pars.target_aperture.upper() == "WFI_CEN":
        tgt_ya_deg, tgt_v2_ref, tgt_v3_ref = get_angles_from_quaternion(q_bam)
    else:
        tgt_sca_info = extract_sca_wcs_params(t_pars.distortion_target_ref)
        tgt_v2_ref = tgt_sca_info["v2_ref"]
        tgt_v3_ref = tgt_sca_info["v3_ref"]
        tgt_ya_deg = tgt_sca_info["v3yangle"]

    tgt_true_ra, tgt_true_dec, true_roll_target, tgt_true_coord = true_point_and_roll(
        tgt_v2_ref, tgt_v3_ref, app_v2v3_frame, velocity
    )
    pa_aperture_true = (true_roll_target + tgt_ya_deg) % 360.0

    # Apparent Reference Coordinates at Target (for DVA anchor)
    tgt_app_coord = SkyCoord(
        lon=tgt_v2_ref * u.arcsec, lat=tgt_v3_ref * u.arcsec, frame=app_v2v3_frame
    ).transform_to("icrs")

    # Velocity Aberration Scale Factor using True Target Vector
    u_tgt_true = tgt_true_coord.cartesian.xyz.value
    _, global_scale_factor = apply_vector_velocity_aberration(u_tgt_true, velocity)

    # 4. SCA-SPECIFIC WCS PROJECTIONS VIA DISTORTION REFERENCE FILE
    true_ra_ref, true_dec_ref, true_roll_ref, _ = true_point_and_roll(
        sca_info["v2_ref"], sca_info["v3_ref"], app_v2v3_frame, velocity
    )

    corners_coord = true_point(sca_info["corners_v2"], sca_info["corners_v3"])

    polygon_str = " ".join(
        [
            f"{ra:.9f} {dec:.9f}"
            for ra, dec in zip(corners_coord.ra.deg, corners_coord.dec.deg, strict=True)
        ]
    )

    # 5. ASDF METADATA ASSEMBLY
    wcs_meta = {
        "pointing": {
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
    obstimes, wcs_metas : [astropy.time.Time[,...]], [dict[...]]
        A 2-tuple is returned with the WCS meta information for the
        the aperture and the V1 axis.
    """
    # Setup structures
    obstimes = []
    wcs_metas = []

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
        wcs_meta = calc_wcs(t_pars)
        obstimes.append(pointing.obstime)
        wcs_metas.append(wcs_meta)

    return obstimes, wcs_metas


def extract_sca_wcs_params(distortion_asdf_path):
    """
    Extracts WCS reference parameters and corner footprints using roman_datamodels.
    Evaluates the linear terms via finite differences (Jacoban) to remain
    agnostic to the ASDF's internal compound model architecture.
    """
    # Active science area is 4088x4088 with reference pixels removed.
    # Corners in 0-indexed space fall exactly on -0.5 and 4087.5
    corners_x_0idx = np.array([-0.5, 4087.5, 4087.5, -0.5])
    corners_y_0idx = np.array([-0.5, -0.5, 4087.5, 4087.5])

    with rdm.open(distortion_asdf_path) as dm:
        detector = str(dm.meta.instrument.detector)
        dist_mod = dm.coordinate_distortion_transform

        # 0-indexed reference pixel for a 4088x4088 array
        x0, y0 = 2043.5, 2043.5

        # 1. Reference coordinates directly from the model
        v2_ref, v3_ref = dist_mod(x0, y0)
        v2_ref, v3_ref = float(v2_ref), float(v3_ref)

        # 2. Finite differences for the Jacobian matrix at the reference pixel
        v2_dx_plus, v3_dx_plus = dist_mod(x0 + 1.0, y0)
        v2_dx_minus, v3_dx_minus = dist_mod(x0 - 1.0, y0)
        dx_v2 = (v2_dx_plus - v2_dx_minus) / 2.0
        dx_v3 = (v3_dx_plus - v3_dx_minus) / 2.0

        v2_dy_plus, v3_dy_plus = dist_mod(x0, y0 + 1.0)
        v2_dy_minus, v3_dy_minus = dist_mod(x0, y0 - 1.0)
        dy_v2 = (v2_dy_plus - v2_dy_minus) / 2.0
        dy_v3 = (v3_dy_plus - v3_dy_minus) / 2.0

        # 3. Derive V3IdlYAngle and VIdlParity from the Jacobian
        v3_idl_yangle = float(np.rad2deg(np.arctan2(dy_v2, dy_v3)))

        # Parity is the sign of the determinant of the [dx, dy] Jacobian matrix
        det_J = (dx_v2 * dy_v3) - (dy_v2 * dx_v3)
        v_idl_parity = int(np.sign(det_J))

        # 4. Evaluate Corner Footprints
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


def v2v3_to_true_icrs(v2_arcsec, v3_arcsec, app_v2v3_frame, velocity):
    """Un-aberrated ICRS positions of V2/V3 points.

    Parameters
    ----------
    v2_arcsec, v3_arcsec : float or array-like
        V2/V3 positions in arcsec.

    app_v2v3_frame : `~astropy.coordinates.SkyOffsetFrame`
        The apparent (aberrated) V2/V3 frame from the attitude quaternion.

    velocity : array-like
        Observatory barycentric velocity (km/s).

    Returns
    -------
    coord : `~astropy.coordinates.SkyCoord`
        True ICRS positions, one per input point.
    """
    app = SkyCoord(
        lon=np.atleast_1d(v2_arcsec) * u.arcsec,
        lat=np.atleast_1d(v3_arcsec) * u.arcsec,
        frame=app_v2v3_frame,
    ).transform_to("icrs")
    vecs = app.cartesian.xyz.value.T
    v_inv = -np.asarray(velocity, dtype=float)
    true = np.array([compute_va_effects_vector(*v_inv, v)[1] for v in vecs])
    ra = np.rad2deg(np.arctan2(true[:, 1], true[:, 0])) % 360.0
    dec = np.rad2deg(np.arcsin(np.clip(true[:, 2], -1.0, 1.0)))
    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")


def true_point_and_roll(v2_arcsec, v3_arcsec, app_v2v3_frame, velocity, step=1.0):
    """True ICRS position of a V2/V3 point and the position angle of +V3 there.

    Parameters
    ----------
    v2_arcsec, v3_arcsec : float
        V2/V3 position in arcsec.

    app_v2v3_frame : `~astropy.coordinates.SkyOffsetFrame`
        The apparent (aberrated) V2/V3 frame from the attitude quaternion.

    velocity : array-like
        Observatory barycentric velocity (km/s).

    step : float
        Step along +V3 (arcsec) used to measure the roll.

    Returns
    -------
    ra, dec, roll, coord : float, float, float, `~astropy.coordinates.SkyCoord`
        True position (deg), +V3 position angle East of North (deg), and the
        position as a scalar SkyCoord.
    """
    c = v2v3_to_true_icrs(
        [v2_arcsec, v2_arcsec],
        [v3_arcsec, v3_arcsec + step],
        app_v2v3_frame,
        velocity,
    )
    roll = float(c[0].position_angle(c[1]).deg % 360.0)
    return float(c[0].ra.deg), float(c[0].dec.deg), roll, c[0]


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
