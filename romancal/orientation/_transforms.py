"""Orientation Direct Cosine Matrix (DCM) transformations and utilities"""

import dataclasses
import logging
import sys
from collections.abc import Callable
from math import cos, sin
from typing import Any

import asdf
import numpy as np
from stcal.velocity_aberration import compute_va_effects_vector

from . import _lib as olib
from . import _pointing as plib

__all__ = []

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LOGLEVELS = [logging.INFO, logging.DEBUG, olib.DEBUG_FULL]

# Define the transformation matrices to move between the Idealized Coordinate System (ICS)
# and the Idealized Coordinate System (Idl). ICS is the spacecraft-centric system used by
# all frames up through the V-frame. Idl is used by the instruments.
# Reference: Innerspace document
M_idl2ics = MX2Z = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
M_ics2idl = MZ2X = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])


# Transforms
@dataclasses.dataclass
class Transforms:
    """The matrices used in calculation of the M_eci2siaf transformation."""

    #: B-frame to FGS
    m_b2fgs: np.ndarray | None = None
    #: ECI to B-frame
    m_eci2b: np.ndarray | None = None
    #: ECI to FGS
    m_eci2fgs: np.ndarray | None = None
    #: ECI to GS
    m_eci2gs: np.ndarray | None = None
    #: ECI to GS apparent
    m_eci2gsapp: np.ndarray | None = None
    #: ECI to V
    m_eci2v: np.ndarray | Any = None
    #: ECI to Guide star apparent position
    m_fgs2gsapp: np.ndarray | None = None
    #: GS apparent to GS corrected
    m_gsapp2gsics: np.ndarray | None = None

    @classmethod
    def from_asdf(cls, asdf_file):
        """
        Create Transforms from AsdfFile.

        Parameters
        ----------
        asdf_file : Stream-like or `asdf.AsdfFile`
            The asdf to create from.

        Returns
        -------
        transforms : Transforms
            The Transforms instance.
        """
        if isinstance(asdf_file, asdf.AsdfFile):
            transforms = asdf_file.tree["transforms"]
        else:
            with asdf.open(asdf_file, memmap=False, lazy_load=False) as af:
                transforms = af.tree["transforms"]

        return cls(**transforms)

    def to_asdf(self):
        """
        Serialize to AsdfFile.

        Returns
        -------
        asdf_file : asdf.AsdfFile
            The ASDF serialization.
        """
        self_dict = dataclasses.asdict(self)
        asdf_file = asdf.AsdfFile({"transforms": self_dict})
        return asdf_file

    def write_to_asdf(self, path):
        """
        Serialize to a file path.

        Parameters
        ----------
        path : Stream-like
            Output file path.
        """
        asdf_file = self.to_asdf()
        asdf_file.write_to(path, all_array_storage="inline")


@dataclasses.dataclass
class TransformParameters:
    """Parameters required for the calculations."""

    #: If telemetry cannot be determined, use existing information in the observation's header.
    allow_default: bool = False
    #: Aperture in use
    aperture: str = ""
    #: Default quaternion to use if engineering is not available.
    default_quaternion: tuple | None = None
    #: Commanded position of the guide star in (H, V) space.
    gscommanded: tuple | None = None
    #: Observation end time
    obsend: float | None = None
    #: Observation start time
    obsstart: float | None = None
    #: The observatory orientation, represented by the ECI quaternion,
    # and other engineering mnemonics
    pointing: plib.Pointing | Any = None
    #: Reduction function to use on values.
    reduce_func: Callable | None = None
    #: Engineering database information
    service_kwargs: dict | None = None
    #: Path to the SIAF folder containing the roman siaf xml definitions
    #: None to use the pysiaf builtin definitions
    siaf_path: str | None = None
    #: If no telemetry can be found during the observation,
    #: the time, in seconds, beyond the observation time to search for telemetry.
    tolerance: float = 60.0
    # Observatory velocity
    velocity: tuple | None = None

    def as_reprdict(self):
        """Return a dict where all values are REPR of their values."""  # numpydoc ignore=RT01
        r = {
            field.name: repr(getattr(self, field.name))
            for field in dataclasses.fields(self)
        }
        return r

    def update_from_engdb(self):
        """Update pointing information."""
        self.pointing = plib.get_pointing(
            self.obsstart,
            self.obsend,
            mnemonics_to_read=plib.COARSE_MNEMONICS,
            service_kwargs=self.service_kwargs,
            tolerance=self.tolerance,
            reduce_func=self.reduce_func,
        )

    def __post_init__(self):
        # Setup the default reduction function.
        if not self.reduce_func:
            self.reduce_func = plib.pointing_from_average


def calc_attitude_matrix(wcs, yangle, position):
    """
    Calculate the DCM attitude from known positions and roll angles.

    Parameters
    ----------
    wcs : WCSRef
        The sky position.

    yangle : float
        The IdlYangle of the point in question.

    position : numpy.array(2)
        The position in Ideal frame.

    Returns
    -------
    m : np.array(3,3)
        The transformation matrix
    """
    # Convert to radians
    ra = wcs.ra * olib.D2R
    dec = wcs.dec * olib.D2R
    yangle_ra = yangle * olib.D2R
    pos_rads = position * olib.A2R
    v2 = pos_rads[0]
    v3 = pos_rads[1]

    # Create the matrices
    r1 = dcm(ra, dec, yangle_ra)

    r2 = np.array(
        [
            [cos(v2) * cos(v3), -sin(v2), -cos(v2) * sin(v3)],
            [sin(v2) * cos(v3), cos(v2), -sin(v2) * sin(v3)],
            [sin(v3), 0.0, cos(v3)],
        ]
    )

    # Final transformation
    m = np.dot(r2, r1)

    logger.debug("attitude DCM: %s", m)
    return m


def calc_gsapp2gs(m_eci2gsapp, velocity):
    """
    Calculate the Velocity Aberration correction.

    This implements Eq. 40 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.2.5:

    The velocity aberration correction is applied in the direction of the guide
    star. The matrix that translates from ECI to the apparent guide star ICS
    frame is M_(ECI→GSAppICS), where the GS Apparent position vector is along
    the z-axis in the guide star ICS frame.

    Note that the convention for Roman is in the opposite direction. The algorithm is the same,
    just different input/output.

    Parameters
    ----------
    m_eci2gsapp : numpy.array(3, 3)
        The the ECI to Guide Star transformation matrix, in the ICS frame.

    velocity : numpy.array([dx, dy, dz])
        The barycentric velocity.

    Returns
    -------
    m_gsapp2gs : numpy.array(3, 3)
        The velocity aberration correction matrix.
    """
    # Check velocity. If present, negate the velocity since
    # the desire is to remove the correction.
    velocity = np.array(velocity)
    if (
        None in velocity
        or np.all(velocity == 0.0)
        or np.sum(np.abs(velocity) > olib.MAX_OBSERVATORY_SPEED)
    ):
        logger.warning(
            "Velocity: %s is either unspecified or contains unreasonable values. Cannot calculate aberration. Returning identity matrix",
            velocity,
        )
        return np.identity(3)
    velocity = -1 * velocity

    # Eq. 35: Guide star position vector
    uz = np.array([0.0, 0.0, 1.0])
    u_gseci = np.dot(np.transpose(m_eci2gsapp), uz)

    # Eq. 36: Compute the apparent shift due to velocity aberration.
    try:
        scale_factor, u_gseci_app = compute_va_effects_vector(*velocity, u_gseci)
    except TypeError:
        logger.warning(
            "Failure in computing velocity aberration. Returning identity matrix."
        )
        logger.warning("Exception: %s", sys.exc_info())
        return np.identity(3)

    # Eq. 39: Rotate from ICS into the guide star frame.
    u_gs_app = np.dot(m_eci2gsapp, u_gseci_app)

    # Eq. 40: Compute the M_gsapp2gs matrix
    u_prod = np.cross(uz, u_gs_app)
    u_prod_mag = np.linalg.norm(u_prod)
    a_hat = u_prod / u_prod_mag
    m_a_hat = np.array(
        [
            [0.0, -a_hat[2], a_hat[1]],
            [a_hat[2], 0.0, -a_hat[0]],
            [-a_hat[1], a_hat[0], 0.0],
        ]
    )
    theta = np.arcsin(u_prod_mag)

    m_gsapp2gs = (
        np.identity(3)
        - (m_a_hat * np.sin(theta))
        + (2 * m_a_hat**2 * np.sin(theta / 2.0) ** 2)
    )

    logger.debug("m_gsapp2gs: %s", m_gsapp2gs)
    return m_gsapp2gs


def calc_m_b2fgs(fgs_q=None):
    """Calculate the B-to-FGS frame DCM

    Parameters
    ----------
    fgs_q : [float, float, float, float] or None
        The quaterion representing the B to FGS transformation.
        If no B-to-FGS quaternion is given, use a pre-launch defined matrix.
    """
    if fgs_q is None:
        logger.warning("No B-to-FGS information is given. Using pre-launch values.")
        fgs_q = olib.FGS_DEFAULT_QUATERNION

    # Calculate the DCM from the quaterion
    return calc_quat2matrix(fgs_q)


def calc_m_fgs2gs(x, y):
    """Calculate DCM that converts FGS reference to guide star location

    Parameters
    ----------
    x, y : float, float
        Guidestar location (radians)
    """
    m = np.array(
        [
            [cos(x), 0.0, sin(x)],
            [-sin(x) * sin(y), cos(y), cos(x) * sin(y)],
            [-sin(x) * cos(y), -sin(y), cos(x) * cos(y)],
        ]
    )

    return m


def calc_m_fgs2gsapp(x, y):
    """Calculate the FGS to Guide star apparent transformation

    Source document: Innerspace Confluence: "Quaternion Transforms for Coarse Pointing WCS"
    from the section "Derive M_fgstogsapp".

    Parameters
    ----------
    x, y : float, float
        Position in arcseconds in the FGS frame

    Returns
    -------
    m_fgs2gsapp : np.array(3,3)
        The DCM representing the transformation.
    """
    x_gs, y_gs = olib.A2R * x, olib.A2R * y
    cx, sx = cos(x_gs), sin(x_gs)
    cy, sy = cos(y_gs), sin(y_gs)

    m_x = np.array([[cx, 0, -sx], [0, 1, 0], [sx, 0, cx]])

    m_y = np.array([[1, 0, 0], [0, cy, sy], [0, -sy, cy]])

    m_gsapp2fgs = np.dot(m_y, m_x)
    m_fgs2gsapp = m_gsapp2fgs.T

    return m_fgs2gsapp


def calc_quat2matrix(q):
    """
    Create a Direction Cosine Matrix (DCM) from a quaternion.

    Parameters
    ----------
    q : np.array(q1, q2, q3, q4)
        Array of quaternions from the engineering database.

    Returns
    -------
    transform : np.array((3, 3))
        The transform matrix representing the transformation
        from observatory orientation to J-Frame.
    """
    q1, q2, q3, q4 = q
    transform = np.array(
        [
            [
                1.0 - 2.0 * q2 * q2 - 2.0 * q3 * q3,
                2.0 * (q1 * q2 + q3 * q4),
                2.0 * (q3 * q1 - q2 * q4),
            ],
            [
                2.0 * (q1 * q2 - q3 * q4),
                1.0 - 2.0 * q3 * q3 - 2.0 * q1 * q1,
                2.0 * (q2 * q3 + q1 * q4),
            ],
            [
                2.0 * (q3 * q1 + q2 * q4),
                2.0 * (q2 * q3 - q1 * q4),
                1.0 - 2.0 * q1 * q1 - 2.0 * q2 * q2,
            ],
        ],
        dtype=float,
    )

    logger.debug("quaternion: %s", transform)
    return transform


def calc_transforms(t_pars: TransformParameters):
    """
    COARSE calculation.

    This implements the overall coarse-mode transfomations as defined by the innerspace
    document "Quaternion Transforms for Coarse Pointing WCS" 2026-01-06.

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms
        The list of coordinate matrix transformations

    Notes
    -----
    """
    logger.info("Calculating transforms...")
    t = Transforms()

    # Quaternion to M_eci2b
    t.m_eci2b = calc_quat2matrix(t_pars.pointing.q)

    # ECI to FGS
    t.m_b2fgs = calc_m_b2fgs(t_pars.pointing.fgs_q)
    t.m_eci2fgs = np.dot(t.m_b2fgs, t.m_eci2b)

    # FGS to Guide star apparent.
    if t_pars.gscommanded is None:
        logger.warning(
            "No command guide star position provided. Assuming guide star is at the aperture reference position."
        )
        hv = (0.0, 0.0)
    else:
        hv = t_pars.gscommanded
    fgs_x, fgs_y = olib.hv_to_fgs(t_pars.aperture, *hv)
    t.m_fgs2gsapp = calc_m_fgs2gsapp(fgs_x, fgs_y)

    # ECI to GS apparent
    t.m_eci2gsapp = np.linalg.multi_dot([t.m_fgs2gsapp, t.m_b2fgs, t.m_eci2b])

    # Use calc_gs2gsapp to convert m_eci2gsapp to VA-applied (or "aberrated") m_eci2gs
    # Note that calc_gs2gsapp should be renamed to calc_gsapp2gs
    t.m_gsapp2gsics = calc_gsapp2gs(t.m_eci2gsapp, t_pars.velocity)

    # ECI to GS
    t.m_eci2gs = np.linalg.multi_dot([M_ics2idl, t.m_gsapp2gsics, t.m_eci2gsapp])

    # ECI to V
    t.m_eci2v = np.linalg.multi_dot(
        [t.m_b2fgs.T, t.m_fgs2gsapp.T, M_idl2ics, t.m_eci2gs]
    )

    return t


def calc_wcs_from_matrix(m):
    """
    Calculate the WCS information from a DCM.

    Parameters
    ----------
    m : np.array((3, 3))
        The DCM matrix to extract WCS information from.

    Returns
    -------
    wcs : WCSRef
        The WCS.
    """
    # V1 RA/Dec is the first row of the transform
    v1_ra, v1_dec = olib.vector_to_angle(m[0])
    wcs = olib.WCSRef(v1_ra, v1_dec, None)

    # V3 is the third row of the transformation
    v3_ra, v3_dec = olib.vector_to_angle(m[2])
    v3wcs = olib.WCSRef(v3_ra, v3_dec, None)

    # Calculate the V3 position angle
    v1_pa = olib.calc_position_angle(wcs, v3wcs)

    # Convert to degrees
    wcs = olib.WCSRef(ra=wcs.ra * olib.R2D, dec=wcs.dec * olib.R2D, pa=v1_pa * olib.R2D)

    logger.debug("wcs: %s", wcs)
    return wcs


def dcm(alpha, delta, angle):
    """
    Construct the Direction Cosine Matrix (DCM).

    Typical usage is passing of (RA, DEC, PositionAngle).
    All values must be in radians.

    Parameters
    ----------
    alpha : float
        First coordinate in radians.

    delta : float
        Second coordinate in radians.

    angle : float
        Position angle in radians.

    Returns
    -------
    dcm : np.array((3, 3))
        The 3x3 direction cosine matrix.
    """
    dcm = np.array(
        [
            [cos(delta) * cos(alpha), cos(delta) * sin(alpha), sin(delta)],
            [
                -cos(angle) * sin(alpha) + sin(angle) * sin(delta) * cos(alpha),
                cos(angle) * cos(alpha) + sin(angle) * sin(delta) * sin(alpha),
                -sin(angle) * cos(delta),
            ],
            [
                -sin(angle) * sin(alpha) - cos(angle) * sin(delta) * cos(alpha),
                sin(angle) * cos(alpha) - cos(angle) * sin(delta) * sin(alpha),
                cos(angle) * cos(delta),
            ],
        ]
    )

    return dcm


def t_pars_from_model(model, t_pars):
    """
    Initialize TransformParameters from a DataModel.

    Parameters
    ----------
    model : DataModel
        Data model to initialize from.

    t_par : TransformParameters
        Transformation parameters updated with model information.
        Updating is performed in-place.
    """
    # Instrument details
    t_pars.aperture = model.meta.wcsinfo.aperture_name
    try:
        exp_type = model.meta.exposure.type.lower()
    except AttributeError:
        exp_type = None
    t_pars.exp_type = exp_type

    # observation parameters
    t_pars.obsstart = model.meta.exposure.start_time
    t_pars.obsend = model.meta.exposure.end_time
    ephem = model.meta.ephemeris
    t_pars.velocity = (ephem.velocity_x, ephem.velocity_y, ephem.velocity_z)

    # Retrieve previously calculated orientation items only if they are currently not defined.
    if t_pars.gscommanded is None:
        t_pars.gscommanded = getattr(model.meta.guide_star, "hv_position", None)
    if t_pars.default_quaternion is None:
        t_pars.default_quaternion = getattr(model.meta.pointing, "quaternion", None)
