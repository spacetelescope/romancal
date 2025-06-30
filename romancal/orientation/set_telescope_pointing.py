"""
Set Telescope Pointing from Observatory Engineering Telemetry.

Calculate and update the pointing-related and world coordinate system-related
keywords. Given a time period, usually defined by an exposure, the engineering
mnemonic database is queried for observatory orientation. The orientation
defines the sky coordinates a particular point on the observatory is pointed to.
Then, using a set of matrix transformations, the sky coordinates of the
reference pixel of a desired aperture is calculated.

The transformations are defined by the Technical Reference JWST-STScI-003222,
SM-12. This document has undergone a number of revisions. The current version
implemented is based on an internal email version Rev. C, produced 2021-11.

There are a number of algorithms, or *methods*, that have been implemented.
Most represent the historical refinement of the algorithm. Until the technical
reference is finalized, all methods will remain in the code. The default,
state-of-the art algorithm is represented by method ``OPS_TR_202111``,
implemented by
`~jwst.lib.set_telescope_pointing.calc_transforms_ops_tr_202111`.

**Interface**

The primary usage is through the command line interface
``set_telescope_pointing.py``. Operating on a list of JWST Level 1b exposures,
this command updates the world coordinate system keywords with the values
necessary to translate from aperture pixel to sky coordinates.

Access to the JWST Engineering Mnemonic database is required. See the
:ref:`Engineering Database Interface <engdb>` for more information.

Programmatically, the command line is implemented by the function
`~jwst.lib.set_telescope_pointing.add_wcs`, which calls the basic function
`~jwst.lib.set_telescope_pointing.calc_wcs`. The available methods are defined
by `~jwst.lib.set_telescope_pointing.Methods`.

There are two data structures used to maintain the state of the transformation.
`~jwst.lib.set_telescope_pointing.TransformParameters` contains the parameters
needed to perform the transformations.
`~jwst.lib.set_telescope_pointing.Transforms` contains the calculated
transformation matrices.

**Transformation Matrices**

All the transformation matrices, as defined by
`~jwst.lib.set_telescope_pointing.Transforms`, are Direction Cosine Matrices
(DCM). A DCM contains the Euler rotation angles that represent the sky
coordinates for a particular frame-of-reference. The initial DCM is provided
through the engineering telemetry and represents where in the sky either the
Fine Guidance Sensor (FGS) or star tracker is pointed to. Then, through a set
of transformations, the DCM for the reference point of the target aperture
is calculated.
"""

import sys

import asdf
from collections import defaultdict, namedtuple
from copy import copy
import dataclasses
from enum import Enum
import logging
from math import cos, sin, sqrt
from typing import Any
from collections.abc import Callable

from astropy.table import Table
from astropy.time import Time
import numpy as np

import roman_datamodels as rdm

from .set_velocity_aberration import compute_va_effects_vector
from ..lib.engdb_tools import ENGDB_Service
from ..lib.siafdb import SIAF, SiafDb

__all__ = [
    "Methods",
    "TransformParameters",
    "Transforms",
    "WCSRef",
    "add_wcs",
    "calc_transforms",
    "calc_wcs",
    "calc_wcs_over_time",
    "update_wcs",
]

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
DEBUG_FULL = logging.DEBUG - 1
LOGLEVELS = [logging.INFO, logging.DEBUG, DEBUG_FULL]

# Datamodels that can be updated, normally
EXPECTED_MODELS = (rdm.datamodels.ScienceRawModel)

# Exposure types that can be updated, normally
TYPES_TO_UPDATE = set()

# Mnemonics for each transformation method.
# dict where value indicates whether the mnemonic is required or not.
COURSE_MNEMONICS_QUATERNION_ECI = [f'SCF_AC_SDR_QBJ_{idx + 1}' for idx in range(4)]
COURSE_MNEMONICS = {q: True for q in COURSE_MNEMONICS_QUATERNION_ECI}

# Conversion from seconds to MJD
SECONDS2MJD = 1 / 24 / 60 / 60

# Conversion of the FCS reference point from the V-Frame.
# This is the pre-launch value, later to be refined and provided
# in the SIAF
M_V2FCS0 = np.array([
    ['-0.0000001', '0.5000141', '0.8660173'],
    ['0.0086567', '-0.8659848', '0.4999953'],
    ['0.9999625', '0.0074969', '-0.0043284']
], dtype=float)

# Default B-frame to FCS frame, M_b_to_fcs
# Pre-launch this is the same as M_v_to_fcs.
M_B2FCS0 = M_V2FCS0

# Define the transformation matrices to move between the Idealized Coordinate System (ICS)
# and the Idealized Coordinate System (Idl). ICS is the spacecraft-centric system used by
# all frames up through the V-frame. Idl is used by the instruments.
# Reference: Eqs. 1 & 2 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
M_idl2ics = MX2Z = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
M_ics2idl = MZ2X = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])

# Degree, radian, angle transformations
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
A2R = D2R / 3600.0
R2A = 3600.0 * R2D
PI2 = np.pi * 2.0


# The available methods for transformation
class Methods(Enum):
    """
    Available methods to calculate V1 and aperture WCS information.
    """
    #: COARSE_NB: Algorithm based solely on the notebook
    COARSE_NB = (
        'coarse_nb',
        'calc_transforms_coarse_nb',
        'calc_wcs_standard',
        COURSE_MNEMONICS
    )

    # Aliases
    #: Algorithm to use by default. Used by Operations.
    default = COARSE_NB
    #: Default algorithm under PCS_MODE COARSE.
    # COARSE = COARSE_TR_202111
    #: Default algorithm for use by Operations.
    OPS = COARSE_NB

    def __new__(cls, value, func_name, calc_func, mnemonics):
        """Create new instance of Methods."""  # numpydoc ignore=RT01
        obj = object.__new__(cls)
        obj._value_ = value
        obj._func_name = func_name  # noqa: SLF001
        obj._calc_func = calc_func  # noqa: SLF001
        obj._mnemonics = mnemonics  # noqa: SLF001
        return obj

    @property
    def calc_func(self):
        """Function associated with the method."""  # numpydoc ignore=RT01
        return globals()[self._calc_func]

    @property
    def func(self):
        """Function associated with the method."""  # numpydoc ignore=RT01
        return globals()[self._func_name]

    @property
    def mnemonics(self):
        """Mnemonics."""  # numpydoc ignore=RT01
        return self._mnemonics

    def __str__(self):
        return self.value


# Pointing container
# Attributes are as follows. Except for the observation time, all values
# are retrieved from the engineering data.
#    obstime      : Time the pointing information refers to.
#    q            : Quaternion of the FGS.
Pointing = namedtuple(
    "Pointing", ['obstime', "q"]
)
Pointing.__new__.__defaults__ = (None,) * 2


# Guide Star ACQ pointing container
# Attributes are as follows. All values are retrieved from the engineering.
#    position : X/Y position of the guide star within the acquisition window of the FGS.
#    corner   : X/Y corner of the acquisition window within the FGS.
#    size     : X/Y size of the acquisition window.
GuideStarPosition = namedtuple("GuideStarPosition", ["position", "corner", "size"])
GuideStarPosition.__new__.__defaults__ = (None,) * 3


# Transforms
@dataclasses.dataclass
class Transforms:
    """The matrices used in calculation of the M_eci2siaf transformation."""

    #: ECI to B-frame
    m_eci2b: np.ndarray | None = None
    #: ECI to FCS
    m_eci2fcs: np.ndarray | None = None
    #: ECI to GS
    m_eci2gs: np.ndarray | None = None
    #: ECI to SIAF
    m_eci2siaf: np.ndarray | None = None
    #: ECI to V
    m_eci2v: np.ndarray | Any = None
    #: Override values. Either another Transforms or dict-like object
    override: object | None = None

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

        Notes
        -----
        The `override` transforms are not serialized, since the values of this transform
        automatically represent what is in the override.
        """
        self_dict = dataclasses.asdict(self)
        del self_dict["override"]  # Do not serialize the override transforms
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

    def __getattribute__(self, name):
        """
        If an override has been specified, return that value regardless.

        Notes
        -----
        This dunder method is called for ALL attributes. Tread carefully.
        """  # numpydoc ignore=RT01
        # If the attribute is not a field, just return its value. Like NOW.
        if name.startswith("_") or name not in self._fields or name == "override":
            return object.__getattribute__(self, name)

        override = self.override
        override_value = getattr(override, name) if override else None
        return override_value if override_value is not None else object.__getattribute__(self, name)

    def __post_init__(self):
        """Post-initialization of a DataClass."""
        # Create a simple list of fields to check against.
        self._fields = [field.name for field in dataclasses.fields(self)]


# WCS reference container
WCSRef = namedtuple("WCSRef", ["ra", "dec", "pa"])
WCSRef.__new__.__defaults__ = (None, None, None)


@dataclasses.dataclass
class TransformParameters:
    """Parameters required for the calculations."""

    #: If telemetry cannot be determined, use existing information in the observation's header.
    allow_default: bool = False
    #: The V3 position angle to use if the pointing information is not found.
    default_pa_v3: float = 0.0
    #: Detector in use.
    detector: str = ""
    #: Do not write out the modified file.
    dry_run: bool = False
    #: URL of the engineering telemetry database REST interface.
    engdb_url: str | None = None
    #: The method, or algorithm, to use in calculating the transform.
    # If not specified, the default method is used.
    method: Methods = Methods.default
    #: Observation end time
    obsend: float | None = None
    #: Observation start time
    obsstart: float | None = None
    #: If set, matrices that should be used instead of the calculated one.
    override_transforms: Transforms | None = None
    #: The observatory orientation, represented by the ECI quaternion,
    # and other engineering mnemonics
    pointing: Pointing | Any = None
    #: Reduction function to use on values.
    reduce_func: Callable | None = None
    #: The SIAF information for the input model
    siaf: SIAF | Any = None
    #: The SIAF database
    siaf_db: SiafDb | Any = None
    #: If no telemetry can be found during the observation,
    #: the time, in seconds, beyond the observation time to search for telemetry.
    tolerance: float = 60.0

    def as_reprdict(self):
        """Return a dict where all values are REPR of their values."""  # numpydoc ignore=RT01
        r = {field.name: repr(getattr(self, field.name)) for field in dataclasses.fields(self)}
        return r

    def update_pointing(self):
        """Update pointing information."""
        self.pointing = get_pointing(
            self.obsstart,
            self.obsend,
            mnemonics_to_read=self.method.mnemonics,
            engdb_url=self.engdb_url,
            tolerance=self.tolerance,
            reduce_func=self.reduce_func,
        )


def add_wcs(
    filename,
    default_pa_v3=0.0,
    siaf_path=None,
    prd=None,
    engdb_url=None,
    tolerance=60,
    allow_default=False,
    reduce_func=None,
    dry_run=False,
    save_transforms=None,
    **transform_kwargs,
):
    """Add WCS information to a JWST DataModel.

    Telescope orientation is attempted to be obtained from
    the engineering database. Failing that, a default pointing
    is used based on proposal target.

    The file is updated in-place.

    Parameters
    ----------
    filename : str
        The path to a data file.

    default_pa_v3 : float
        The V3 position angle to use if the pointing information
        is not found.

    siaf_path : str or file-like object or None
        The path to the SIAF database. See `SiafDb` for more information.

    prd : str
        The PRD version from the `pysiaf` to use.
        `siaf_path` overrides this value.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default : bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func : func or None
        Reduction function to use on values.

    dry_run : bool
        Do not write out the modified file.

    save_transforms : Path-like or None
        File to save the calculated transforms to.

    **transform_kwargs : dict
        Keyword arguments used by matrix calculation routines.

    Notes
    -----
    This function adds absolute pointing information to the Roman datamodels
    provided. By default, only Stage 1 exposures are allowed to be updated.
    These have the suffixes of "uncal" representing datamodels ScienceRawModel.
    Any higher level product, from Stage 2b and beyond, that has had the
    `assign_wcs` step applied, have improved WCS information. Running this task
    on such files will potentially corrupt the WCS.
    """
    logger.info("Updating WCS info for file %s", filename)

    with rdm.open(filename) as model:
        if not isinstance(model, EXPECTED_MODELS):
            logger.warning("Input %s is not of an expected type (uncal)", model)
            logger.warning(
                "    Updating pointing may have no effect or detrimental effects on the "
                "WCS information,"
            )
            logger.warning(
                "    especially if the input is the result of Level2b or higher calibration."
            )
            raise TypeError(
                f"Input model {model} is not one of {EXPECTED_MODELS}."
                "\n\tFailing WCS processing."
            )

        t_pars, transforms = update_wcs(
            model,
            default_pa_v3=default_pa_v3,
            siaf_path=siaf_path,
            prd=prd,
            engdb_url=engdb_url,
            tolerance=tolerance,
            allow_default=allow_default,
            reduce_func=reduce_func,
            **transform_kwargs,
        )

        if dry_run:
            logger.info("Dry run requested; results are not saved.")
        else:
            logger.info("Saving updated model %s", filename)
            model.save(filename)
            if transforms and save_transforms:
                logger.info("Saving transform matrices to %s", save_transforms)
                transforms.write_to_asdf(save_transforms)

    logger.info("...update completed")

def update_wcs(
    model,
    default_pa_v3=0.0,
    default_roll_ref=0.0,
    siaf_path=None,
    prd=None,
    engdb_url=None,
    tolerance=60,
    allow_default=False,
    reduce_func=None,
    **transform_kwargs,
):
    """
    Update WCS pointing information.

    Given a `jwst.datamodels.JwstDataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        The model to update.

    default_roll_ref : float
        If pointing information cannot be retrieved,
        use this as the roll ref angle.

    siaf_path : str or Path-like object
        The path to the SIAF database. See `SiafDb` for more information.

    prd : str
        The PRD version from the `pysiaf` to use.
        `siaf_path` overrides this value.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    allow_default : bool
        If telemetry cannot be determine, use existing
        information in the observation's header.

    reduce_func : func or None
        Reduction function to use on values.

    **transform_kwargs : dict
        Keyword arguments used by matrix calculation routines.

    Returns
    -------
    t_pars, transforms : TransformParameters, Transforms
        The parameters and transforms calculated. May be
        None for either if telemetry calculations were not
        performed. In particular, FGS GUIDER data does
        not need `transforms`.
    """
    t_pars = transforms = None  # Assume telemetry is not used.

    # TODO: determine how to find prd versions
    # if not prd:
    #     prd = model.meta.prd_version
    # siaf_db = SiafDb(source=siaf_path, prd=prd)

    # Configure transformation parameters.
    t_pars = t_pars_from_model(
        model,
        default_pa_v3=default_pa_v3,
        engdb_url=engdb_url,
        tolerance=tolerance,
        allow_default=allow_default,
        reduce_func=reduce_func,
        **transform_kwargs,
    )

    # TODO: Determine if this is necessary
    # Populate header with SIAF information.
    # populate_model_from_siaf(model, t_pars.siaf)

    # Calculate WCS.
    transforms = update_wcs_from_telem(model, t_pars)

    return t_pars, transforms


def update_wcs_from_telem(model, t_pars: TransformParameters):
    """
    Update WCS pointing information.

    Given a `jwst.datamodels.JwstDataModel`, determine the simple WCS parameters
    from the SIAF keywords in the model and the engineering parameters
    that contain information about the telescope pointing.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    model : `~jwst.datamodels.JwstDataModel`
        The model to update. The update is done in-place.

    t_pars : `TransformParameters`
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : Transforms or None
        If available, the transformation matrices.
    """
    logger.info("Updating wcs from telemetry.")
    transforms = None  # Assume no transforms are calculated.

    # Setup default WCS info if actual pointing and calculations fail.
    wcsinfo = WCSRef(model.meta.pointing.target_ra, model.meta.pointing.target_dec, model.meta.pointing.pa_v3)
    vinfo = wcsinfo

    # Get the pointing information
    try:
        t_pars.update_pointing()
    except ValueError as exception:
        if not t_pars.allow_default:
            raise
        else:
            logger.warning(
                "Cannot retrieve valid telescope pointing."
                " Default pointing parameters will be used."
            )
            logger.warning("Exception is %s", exception)
            logger.info("Setting ENGQLPTG keyword to PLANNED")
            model.meta.visit.engdb_pointing_quality = "PLANNED"
    else:
        logger.info("Successful read of engineering quaternions:")
        logger.info("\tPointing: %s", t_pars.pointing)

    # If pointing is available, attempt to calculate WCS information
    if t_pars.pointing is not None:
        try:
            wcsinfo, vinfo, transforms = calc_wcs(t_pars)
            pointing_engdb_quality = f"CALCULATED_{t_pars.method.value.upper()}"
            logger.info("Setting ENGQLPTG keyword to %s", pointing_engdb_quality)
            model.meta.visit.engdb_pointing_quality = pointing_engdb_quality
        except Exception as e:
            logger.warning(
                "WCS calculation has failed and will be skipped."
                "Default pointing parameters will be used."
            )
            logger.warning("Exception is %s", e)
            if not t_pars.allow_default:
                raise
            else:
                logger.info("Setting ENGQLPTG keyword to PLANNED")
                model.meta.visit.engdb_pointing_quality = "PLANNED"
    logger.info("Aperture WCS info: %s", wcsinfo)
    logger.info("V1 WCS info: %s", vinfo)

    # Update V1 pointing
    model.meta.pointing.ra_v1 = vinfo.ra
    model.meta.pointing.dec_v1 = vinfo.dec
    model.meta.pointing.pa_v3 = vinfo.pa

    # Update Aperture pointing
    model.meta.pointing.pa_aperture = wcsinfo.pa
    model.meta.wcsinfo.ra_ref = wcsinfo.ra
    model.meta.wcsinfo.dec_ref = wcsinfo.dec

    # TODO: when siaf info is truly incorporated
    # model.meta.wcsinfo.roll_ref = pa_to_roll_ref(wcsinfo.pa, t_pars.siaf)
    model.meta.wcsinfo.roll_ref = wcsinfo.pa

    # TODO: when siaf info is truly incorporated
    # Calculate S_REGION with the footprint
    # information
    # try:
    #     update_s_region(model, t_pars.siaf)
    # except Exception as e:
    #     logger.warning("Calculation of S_REGION failed and will be skipped.")
    #     logger.warning("Exception is %s", e)

    return transforms


def calc_wcs_over_time(obsstart, obsend, t_pars: TransformParameters):
    """
    Calculate V1 and WCS over a time period.

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
        pointings = get_pointing(
            obsstart,
            obsend,
            engdb_url=t_pars.engdb_url,
            tolerance=t_pars.tolerance,
            reduce_func=t_pars.reduce_func,
        )
    except ValueError:
        logger.warning("Cannot get valid engineering mnemonics from engineering database")
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


def calc_wcs(t_pars: TransformParameters):
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
    if t_pars.siaf is None:
        t_pars.siaf = SIAF()

    # Calculate transforms
    transforms = calc_transforms(t_pars)

    # Calculate the wcs information
    wcsinfo, vinfo = t_pars.method.calc_func(transforms)

    # That's all folks
    return wcsinfo, vinfo, transforms


def calc_wcs_standard(transforms: Transforms):
    """
    Calculate WCS transformation.

    Given observatory orientation and target aperture,
    calculate V1 and Reference Pixel sky coordinates.

    Parameters
    ----------
    transforms : Transforms
        The transformation matrices.

    Returns
    -------
    wcsinfo, vinfo : WCSRef, WCSRef
        A 2-tuple is returned with the WCS pointing for
        the aperture and the V1 axis.
    """
    # Calculate the V1 WCS information
    vinfo = calc_wcs_from_matrix(transforms.m_eci2v)

    # Calculate the Aperture WCS
    wcsinfo = calc_wcs_from_matrix(transforms.m_eci2siaf)

    # That's all folks
    return wcsinfo, vinfo


def calc_transforms(t_pars: TransformParameters):
    """
    Calculate transforms  which determine reference point celestial WCS.

    This implements Eq. 3 from Technical Report JWST-STScI-003222, SM-12. Rev. C, 2021-11
    From Section 3:

    The Direction Cosine Matrix (DCM) that provides the transformation of a
    unit pointing vector defined in inertial frame (ECI J2000) coordinates to a
    unit vector defined in the science aperture Ideal frame coordinates is
    defined as [follows.]

    Parameters
    ----------
    t_pars : TransformParameters
        The transformation parameters. Parameters are updated during processing.

    Returns
    -------
    transforms : `Transforms`
        The list of coordinate matrix transformations.
    """
    t_pars.method = t_pars.method if t_pars.method else Methods.default

    transforms = t_pars.method.func(t_pars)
    return transforms


def calc_transforms_coarse_nb(t_pars: TransformParameters):
    """
    COARSE calculation.

    This implements the overall coarse-mode transfomations as defined by the
    example notebook as provided by T.Sohn.

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
    logger.info(
        "Calculating transforms using COARSE_NB Tracking..."
    )
    t_pars.method = Methods.COARSE_NB
    t = Transforms(override=t_pars.override_transforms)

    # Quaternion to M_eci2b
    t.m_eci2b = calc_m_eci2b(t_pars.pointing.q)

    # ECI to FCS
    t.m_eci2fcs = np.dot(M_B2FCS0, t.m_eci2b)

    # ECI to GS
    t.m_eci2gs = np.dot(M_ics2idl, t.m_eci2fcs)

    # ECI to V
    t.m_eci2v = np.linalg.multi_dot([M_B2FCS0.T, M_idl2ics, t.m_eci2gs])

    # ECI to SIAF
    t.m_eci2siaf = np.linalg.multi_dot([M_ics2idl, M_B2FCS0, t.m_eci2v])

    return t

def calc_gs2gsapp(m_eci2gsics, jwst_velocity):
    """
    Calculate the Velocity Aberration correction.

    This implements Eq. 40 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.2.5:

    The velocity aberration correction is applied in the direction of the guide
    star. The matrix that translates from ECI to the apparent guide star ICS
    frame is M_(ECI→GSAppICS), where the GS Apparent position vector is along
    the z-axis in the guide star ICS frame.

    Parameters
    ----------
    m_eci2gsics : numpy.array(3, 3)
        The the ECI to Guide Star transformation matrix, in the ICS frame.

    jwst_velocity : numpy.array([dx, dy, dz])
        The barycentric velocity of JWST.

    Returns
    -------
    m_gs2gsapp : numpy.array(3, 3)
        The velocity aberration correction matrix.
    """
    # Check velocity. If present, negate the velocity since
    # the desire is to remove the correction.
    if jwst_velocity is None or any(jwst_velocity == None):  # noqa: E711 Syntax needed for numpy arrays.
        logger.warning(
            "Velocity: %s contains None. Cannot calculate aberration. Returning identity matrix",
            jwst_velocity,
        )
        return np.identity(3)
    velocity = -1 * jwst_velocity

    # Eq. 35: Guide star position vector
    uz = np.array([0.0, 0.0, 1.0])
    u_gseci = np.dot(np.transpose(m_eci2gsics), uz)

    # Eq. 36: Compute the apparent shift due to velocity aberration.
    try:
        scale_factor, u_gseci_app = compute_va_effects_vector(*velocity, u_gseci)
    except TypeError:
        logger.warning("Failure in computing velocity aberration. Returning identity matrix.")
        logger.warning("Exception: %s", sys.exc_info())
        return np.identity(3)

    # Eq. 39: Rotate from ICS into the guide star frame.
    u_gs_app = np.dot(m_eci2gsics, u_gseci_app)

    # Eq. 40: Compute the M_gs2gsapp matrix
    u_prod = np.cross(uz, u_gs_app)
    u_prod_mag = np.linalg.norm(u_prod)
    a_hat = u_prod / u_prod_mag
    m_a_hat = np.array(
        [[0.0, -a_hat[2], a_hat[1]], [a_hat[2], 0.0, -a_hat[0]], [-a_hat[1], a_hat[0], 0.0]]
    )
    theta = np.arcsin(u_prod_mag)

    m_gs2gsapp = (
        np.identity(3) - (m_a_hat * np.sin(theta)) + (2 * m_a_hat**2 * np.sin(theta / 2.0) ** 2)
    )

    logger.debug("m_gs2gsapp: %s", m_gs2gsapp)
    return m_gs2gsapp


def calc_attitude_matrix(wcs, yangle, position):
    """
    Calculate the DCM attitude from known positions and roll angles.

    This implements Appendix A from Technical Report JWST-STScI-003222, SM-12. 2021-07.

    Parameters
    ----------
    wcs : WCSRef
        The guide star position.

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
    ra = wcs.ra * D2R
    dec = wcs.dec * D2R
    yangle_ra = yangle * D2R
    pos_rads = position * A2R
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
    v1_ra, v1_dec = vector_to_angle(m[0])
    wcs = WCSRef(v1_ra, v1_dec, None)

    # V3 is the third row of the transformation
    v3_ra, v3_dec = vector_to_angle(m[2])
    v3wcs = WCSRef(v3_ra, v3_dec, None)

    # Calculate the V3 position angle
    v1_pa = calc_position_angle(wcs, v3wcs)

    # Convert to degrees
    wcs = WCSRef(ra=wcs.ra * R2D, dec=wcs.dec * R2D, pa=v1_pa * R2D)

    logger.debug("wcs: %s", wcs)
    return wcs


def calc_m_eci2b(q):
    """
    Calculate ECI to B-frame matrix from quaternions.


    This implements the M_eci_to_b calculation as presented in the
    STScI Innerspace document "Quaternion Transforms for Coarse Pointing WCS".

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
        dtype=float
    )

    logger.debug("quaternion: %s", transform)
    return transform


def calc_v2siaf_matrix(siaf):
    """
    Calculate the SIAF transformation matrix.

    This implements Eq. 12 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.1:

    The V to SIAF parameters V3IdlYang, V2Ref, V3Ref, and VIdlParity are
    defined and their usage explained in SIAF2017. The parameter values for
    each aperture are specified in the Project Reference Database (PRD).

    Parameters
    ----------
    siaf : SIAF
        The SIAF parameters, where angles are in arcseconds/degrees.

    Returns
    -------
    transform : np.array((3, 3))
        The V1 to SIAF transformation matrix.
    """
    v2, v3, v3idlyang, vparity = (siaf.v2_ref, siaf.v3_ref, siaf.v3yangle, siaf.vparity)
    mat = dcm(v2 * A2R, v3 * A2R, v3idlyang * D2R)
    pmat = np.array([[0.0, vparity, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])

    transform = np.dot(pmat, mat)

    logger.debug("transform: %s", transform)
    return transform


def calc_position_angle(point, ref):
    """
    Calculate position angle from reference to point.

    Algorithm implemented is from JWST Technical Report JWST-STScI-001550, SM-12,
    2017-11-08, Rev A., Section 5.2, page 29, final equation::

        tan(pa) = cos(dec_r) * sin(ra_r - ra_p) / (sin(dec_r)cos(dec_p) - cos(dec_r)sin(dec_p)cos(ra_r-ra_p))

    where::

        pa : position angle
        *_r : reference
        *_p : point

    Typically the reference is the V3 RA/DEC and point is the object RA/DEC.

    Parameters
    ----------
    point : WCSRef
        The POINT wcs parameters, in radians.

    ref : WCSRef
        The TARGET wcs parameters, in radians.

    Returns
    -------
    point_pa : float
      The POINT position angle, in radians
    """  # noqa: E501
    y = cos(ref.dec) * sin(ref.ra - point.ra)
    x = sin(ref.dec) * cos(point.dec) - cos(ref.dec) * sin(point.dec) * cos(ref.ra - point.ra)
    point_pa = np.arctan2(y, x)
    if point_pa < 0:
        point_pa += PI2
    if point_pa >= PI2:
        point_pa -= PI2

    logger.debug("Given reference: %s, point: %s, then PA: %s", ref, point, point_pa)
    return point_pa


def get_pointing(
        obsstart,
        obsend,
        mnemonics_to_read=COURSE_MNEMONICS,
        engdb_url=None,
        tolerance=60,
        reduce_func=None,
):
    """
    Get telescope pointing engineering data.

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    mnemonics_to_read : {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    reduce_func : func or None
        Reduction function to use on values.
        If None, the average pointing is returned.

    Returns
    -------
    pointing : Pointing or [Pointing(, ...)]
        The engineering pointing parameters.
        If the `result_type` is `all`, a list
        of pointings will be returned.

    Raises
    ------
    ValueError
        Cannot retrieve engineering information.

    Notes
    -----
    For the moment, the first found values will be used.
    This will need be re-examined when more information is
    available.
    """
    if reduce_func is None:
        reduce_func = pointing_from_average

    logger.info("Determining pointing between observations times (mjd):")
    logger.info("obsstart: %s obsend: %s", obsstart, obsend)
    logger.info("Telemetry search tolerance: %s", tolerance)
    logger.info("Reduction function: %s", reduce_func)

    mnemonics = get_mnemonics(
        obsstart,
        obsend,
        mnemonics_to_read=mnemonics_to_read,
        tolerance=tolerance,
        engdb_url=engdb_url,
    )
    reduced = reduce_func(mnemonics_to_read, mnemonics)

    logger.log(DEBUG_FULL, "Mnemonics found:")
    logger.log(DEBUG_FULL, "%s", mnemonics)
    logger.info("Reduced set of pointings:")
    logger.info("%s", reduced)

    return reduced


def vector_to_angle(v):
    """
    Return tuple of spherical angles from unit direction Vector.

    This implements Eq. 10 & 11 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3:

    The Direction Cosine Matrix (DCM) that provides the transformation of a
    unit pointing vector defined in inertial frame (ECI J2000) coordinates to a
    unit vector defined in the science aperture Ideal frame coordinates is
    defined as.

    Parameters
    ----------
    v : [v0, v1, v2]
        Direction vector.

    Returns
    -------
    alpha, delta : float, float
        The spherical angles, in radians.
    """
    alpha = np.arctan2(v[1], v[0])
    delta = np.arcsin(v[2])
    if alpha < 0.0:
        alpha += 2.0 * np.pi
    return alpha, delta


def angle_to_vector(alpha, delta):
    """
    Convert spherical angles to unit vector.

    This implements Eq. 9 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3.

    Parameters
    ----------
    alpha, delta : float
        Spherical angles in radians.

    Returns
    -------
    v : [float, float, float]
        Unit vector.
    """
    v0 = cos(delta) * cos(alpha)
    v1 = cos(delta) * sin(alpha)
    v2 = sin(delta)

    return [v0, v1, v2]


def compute_local_roll(pa_v3, ra_ref, dec_ref, v2_ref, v3_ref):
    """
    Compute local roll.

    Compute the position angle of V3 (measured N to E)
    at the center af an aperture.

    Parameters
    ----------
    pa_v3 : float
        Position angle of V3 at (V2, V3) = (0, 0) [in deg].
    ra_ref, dec_ref : float
        RA and DEC corresponding to V2_REF and V3_REF, [in deg].
    v2_ref, v3_ref : float
        Reference point in the V2, V3 frame [in arcsec].

    Returns
    -------
    new_roll : float
        The value of ROLL_REF (in deg).
    """
    v2 = np.deg2rad(v2_ref / 3600)
    v3 = np.deg2rad(v3_ref / 3600)
    ra_ref = np.deg2rad(ra_ref)
    dec_ref = np.deg2rad(dec_ref)
    pa_v3 = np.deg2rad(pa_v3)

    m = np.array(
        [
            [
                cos(ra_ref) * cos(dec_ref),
                -sin(ra_ref) * cos(pa_v3) + cos(ra_ref) * sin(dec_ref) * sin(pa_v3),
                -sin(ra_ref) * sin(pa_v3) - cos(ra_ref) * sin(dec_ref) * cos(pa_v3),
            ],
            [
                sin(ra_ref) * cos(dec_ref),
                cos(ra_ref) * cos(pa_v3) + sin(ra_ref) * sin(dec_ref) * sin(pa_v3),
                cos(ra_ref) * sin(pa_v3) - sin(ra_ref) * sin(dec_ref) * cos(pa_v3),
            ],
            [sin(dec_ref), -cos(dec_ref) * sin(pa_v3), cos(dec_ref) * cos(pa_v3)],
        ]
    )

    return _roll_angle_from_matrix(m, v2, v3)


def _roll_angle_from_matrix(matrix, v2, v3):
    x = -(matrix[2, 0] * np.cos(v2) + matrix[2, 1] * np.sin(v2)) * np.sin(v3) + matrix[
        2, 2
    ] * np.cos(v3)
    y = (matrix[0, 0] * matrix[1, 2] - matrix[1, 0] * matrix[0, 2]) * np.cos(v2) + (
        matrix[0, 1] * matrix[1, 2] - matrix[1, 1] * matrix[0, 2]
    ) * np.sin(v2)
    new_roll = np.rad2deg(np.arctan2(y, x))
    if new_roll < 0:
        new_roll += 360
    return new_roll


def get_mnemonics(
    obsstart, obsend, tolerance, mnemonics_to_read=COURSE_MNEMONICS, engdb_url=None
):
    """
    Retrieve pointing mnemonics from the engineering database.

    Parameters
    ----------
    obsstart, obsend : float
        astropy.Time observation start/end times.

    tolerance : int
        If no telemetry can be found during the observation,
        the time, in seconds, beyond the observation time to
        search for telemetry.

    mnemonics_to_read : {str: bool[,...]}
        The mnemonics to fetch. key is the mnemonic and
        value is whether it is required to be found.

    engdb_url : str or None
        URL of the engineering telemetry database REST interface.

    Returns
    -------
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Raises
    ------
    ValueError
        Cannot retrieve engineering information.
    """
    try:
        engdb = ENGDB_Service(base_url=engdb_url)
    except Exception as exception:
        raise ValueError(f"Cannot open engineering DB connection\nException: {exception}") from None
    logger.info("Querying engineering DB: %s", engdb.base_url)

    # Construct the mnemonic values structure.
    mnemonics = {mnemonic: None for mnemonic in mnemonics_to_read}

    # Retrieve the mnemonics from the engineering database.
    # Check for whether the bracket values are used and
    # within tolerance.
    for mnemonic in mnemonics:
        try:
            mnemonics[mnemonic] = engdb.get_values(
                mnemonic,
                obsstart,
                obsend,
                time_format="mjd",
                include_obstime=True,
                include_bracket_values=False,
            )
        except Exception as exception:
            raise ValueError(f"Cannot retrieve {mnemonic} from engineering.") from exception

        # If more than two points exist, throw off the bracket values.
        # Else, ensure the bracket values are within the allowed time.
        if len(mnemonics[mnemonic]) < 2:
            logger.warning("Mnemonic %s has no telemetry within the observation time.", mnemonic)
            logger.warning("Attempting to use bracket values within %s seconds", tolerance)

            mnemonics[mnemonic] = engdb.get_values(
                mnemonic,
                obsstart,
                obsend,
                time_format="mjd",
                include_obstime=True,
                include_bracket_values=True,
            )

            tolerance_mjd = tolerance * SECONDS2MJD
            allowed_start = obsstart - tolerance_mjd
            allowed_end = obsend + tolerance_mjd
            allowed = [
                value
                for value in mnemonics[mnemonic]
                if allowed_start <= value.obstime <= allowed_end
            ]
            if not len(allowed):
                raise ValueError(
                    "No telemetry exists for mnemonic {} within {} and {}".format(
                        mnemonic,
                        Time(allowed_start, format="mjd").isot,
                        Time(allowed_end, format="mjd").isot,
                    )
                )
            mnemonics[mnemonic] = allowed

    # All mnemonics must have some values.
    if not all(len(mnemonic) for mnemonic in mnemonics.values()):
        raise ValueError("Incomplete set of pointing mnemonics")

    return mnemonics


def all_pointings(mnemonics_to_read, mnemonics):  # noqa: ARG001
    """
    V1 of making pointings.

    Parameters
    ----------
    mnemonics_to_read : {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Returns
    -------
    pointings : [Pointing[,...]]
        List of pointings.
    """
    pointings = []
    filled = fill_mnemonics_chronologically(mnemonics)
    for obstime, mnemonics_at_time in filled.items():
        # Fill out the matrices
        q = np.array([mnemonics_at_time[m].value for m in COURSE_MNEMONICS_QUATERNION_ECI])

        pointing = Pointing(
            q=q,
            obstime=obstime,
        )
        pointings.append(pointing)

    if not len(pointings):
        raise ValueError("No non-zero quaternion found.")

    return pointings


def first_pointing(mnemonics_to_read, mnemonics):
    """
    Return first pointing.

    Parameters
    ----------
    mnemonics_to_read : {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Returns
    -------
    pointing : Pointing
        First pointing.
    """
    pointings = all_pointings(mnemonics_to_read, mnemonics)
    return pointings[0]


def pointing_from_average(mnemonics_to_read, mnemonics):
    """
    Determine single pointing from average of available pointings.

    Parameters
    ----------
    mnemonics_to_read : {str: bool[,...]}
        The mnemonics to read. Key is the mnemonic name.
        Value is a boolean indicating whether the mnemonic
        is required to have values or not.

    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Returns
    -------
    pointing : Pointing
        Pointing from average.
    """
    # Get average observation time.
    times = [
        eng_param.obstime.unix
        for key in mnemonics
        for eng_param in mnemonics[key]
        if eng_param.obstime.unix != 0.0
    ]
    if len(times) > 0:
        obstime = Time(np.average(times), format="unix")
    else:
        raise ValueError("No valid times in range")

    # Get averages for all the mnemonics.
    mnemonic_averages = {}
    zero_mnemonics = []
    for mnemonic in mnemonics:
        values = [eng_param.value for eng_param in mnemonics[mnemonic]]
        # Weed out mnemonic entries that are zero, though some are OK to be zero.
        if mnemonics_to_read[mnemonic]:
            good_mnemonic = []
            for this_value in values:
                if this_value != 0.0:
                    good_mnemonic.append(this_value)
            if len(good_mnemonic) > 0:
                mnemonic_averages[mnemonic] = np.average(good_mnemonic)
            else:
                zero_mnemonics.append(mnemonic)
        else:
            mnemonic_averages[mnemonic] = np.average(values)

    # Raise exception if there are mnemonics with only zeros in the time range
    if len(zero_mnemonics):
        logger.warning(
            "The following engineering mnemonics only contained zeros in the requested "
            "time interval:"
        )
        badmnemonicsstring = " ".join(zero_mnemonics)
        logger.info(badmnemonicsstring)
        raise ValueError("Bad telemetry values")

    # Fill out the pointing matrices.
    q = np.array([mnemonic_averages[m] for m in COURSE_MNEMONICS_QUATERNION_ECI])

    pointing = Pointing(
        obstime=obstime,
        q=q,
    )

    # That's all folks
    return pointing


def fill_mnemonics_chronologically(mnemonics, filled_only=True):
    """
    Return time-ordered mnemonic list with progressive values.

    The different set of mnemonics used for observatory orientation
    appear at different cadences. This routine creates a time-ordered dictionary
    with all the mnemonics for each time found in the engineering. For mnemonics
    missing for a particular time, the last previous value is used.

    Parameters
    ----------
    mnemonics : {mnemonic: [value[,...]]}
        Dictionary mapping mnemonics to their respective values.

    filled_only : bool
        Only return a matrix where observation times have all the mnemonics defined.

    Returns
    -------
    filled_by_time : {obstime: {mnemonic: value}}
        Time-ordered mnemonic list with progressive values.
    """
    # Collect all information by observation time and order.
    by_obstime = defaultdict(dict)
    n_mnemonics = len(mnemonics)
    for mnemonic, values in mnemonics.items():
        for value in values:
            by_obstime[value.obstime][mnemonic] = value
    by_obstime = sorted(by_obstime.items())

    # Created the filled matrix
    filled = {}
    last_obstime = {}
    for obstime, mnemonics_at_time in by_obstime:
        last_obstime.update(mnemonics_at_time)
        if len(last_obstime) >= n_mnemonics or not filled_only:
            # Engineering data may be present, but all zeros.
            # Filter out this situation also.
            if filled_only:
                values = [value.value for value in last_obstime.values()]
                if not any(values):
                    continue

            filled[obstime] = copy(last_obstime)

    return filled


def fill_mnemonics_chronologically_table(mnemonics, filled_only=True):
    """
    Return time-ordered mnemonic list with progressive values.

    The different set of mnemonics used for observatory orientation
    appear at different cadences. This routine creates a time-ordered dictionary
    with all the mnemonics for each time found in the engineering. For mnemonics
    missing for a particular time, the last previous value is used.

    Parameters
    ----------
    mnemonics : {mnemonic: [value[,...]]}
        Dictionary mapping mnemonics to their respective values.

    filled_only : bool
        Only return a matrix where observation times have all the mnemonics defined.

    Returns
    -------
    filled_by_time : `astropy.table.Table`
        Time-ordered mnemonic list with progressive values.
    """
    filled = fill_mnemonics_chronologically(mnemonics, filled_only=filled_only)

    names = list(mnemonics.keys())
    names = ["time"] + names
    time_idx = 0

    values = [[] for _ in names]

    for time in filled:
        values[time_idx].append(time)
        for mnemonic in filled[time]:
            idx = names.index(mnemonic)
            values[idx].append(filled[time][mnemonic].value)

    t = Table(values, names=names)

    return t


def position_to_dcm(x, y):
    """
    Calculate the Direction Cosine Matrix for a given X,Y position.

    Parameters
    ----------
    x, y : float, float
        Position in arcseconds

    Returns
    -------
    dcm : np.array(size=(3, 3))
        The direction cosine matrix
    """
    dcm = np.array([
        [cos(x),           0,       sin(x)],
        [-sin(x) * sin(y), cos(y),  cos(x) * sin(y)],
        [-sin(x) * cos(y), -sin(y), cos(x) * cos(y)]
    ])

    return dcm


def cart_to_vector(coord):
    """
    Convert Cartesian to a unit vector.

    This implements Eq. 6 from Technical Report JWST-STScI-003222, SM-12, Rev. C, 2021-11
    From Section 3:

    The Direction Cosine Matrix (DCM) that provides the transformation of a
    unit pointing vector defined in inertial frame (ECI J2000) coordinates to a
    unit vector defined in the science aperture Ideal frame coordinates is
    defined as.

    Parameters
    ----------
    coord : numpy.array(2)
        The Cartesian coordinate.

    Returns
    -------
    vector : numpy.array(3)
        The vector version.
    """
    vector = np.array([coord[0], coord[1], sqrt(1 - coord[0] ** 2 - coord[1] ** 2)])

    return vector


def pa_to_roll_ref(pa: float, siaf: SIAF):
    """
    Calculate Roll from the position angle of the given aperture.

    Parameters
    ----------
    pa : float
        Position angle of the aperture, in degrees.

    siaf : SIAF
        The SIAF of the aperture.

    Returns
    -------
    roll_ref : float
        The roll reference, in degrees.
    """
    return pa - siaf.v3yangle


def t_pars_from_model(model, **t_pars_kwargs):
    """
    Initialize TransformParameters from a DataModel.

    Parameters
    ----------
    model : DataModel
        Data model to initialize from.

    **t_pars_kwargs : dict
        Keyword arguments used to initialize the TransformParameters object
        before reading from the model meta information.

    Returns
    -------
    t_par : TransformParameters
        The initialized parameters.
    """
    t_pars = TransformParameters(**t_pars_kwargs)

    # Retrieve SIAF information
    if t_pars.siaf is None:
        siaf = None
        useafter = None
        if t_pars.siaf_db is not None:
            aperture_name = model.meta.aperture.name.upper()
            useafter = model.meta.observation.date
            if aperture_name != "UNKNOWN":
                logger.info("Updating WCS for aperture %s", aperture_name)

                # Special case. With aperture MIRIM_TAMRS, the siaf definition is
                # for the subarray of interest. However, the whole detector is
                # read out. Hence, need to convert pixel coordinates to be detector-based.
                to_detector = False
                if aperture_name == "MIRIM_TAMRS":
                    to_detector = True
                siaf = t_pars.siaf_db.get_wcs(
                    aperture_name, to_detector=to_detector, useafter=useafter
                )
        t_pars.siaf = siaf
        t_pars.useafter = useafter
    logger.debug("SIAF: %s", t_pars.siaf)

    # Instrument details
    t_pars.detector = model.meta.instrument.detector
    try:
        exp_type = model.meta.exposure.type.lower()
    except AttributeError:
        exp_type = None
    t_pars.exp_type = exp_type

    # observation parameters
    t_pars.obsstart = model.meta.exposure.start_time
    t_pars.obsend = model.meta.exposure.end_time
    logger.debug("Observation time: %s - %s", t_pars.obsstart, t_pars.obsend)

    # Set the transform and WCS calculation method.
    t_pars.method = Methods.default

    # Set pointing reduction function if not already set.
    if not t_pars.reduce_func:
        t_pars.reduce_func = pointing_from_average

    return t_pars


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
