"""Transformation parameters used in the orientation calculations."""

import dataclasses
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from . import _lib as olib
from . import _pointing as plib

__all__ = []

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LOGLEVELS = [logging.INFO, logging.DEBUG, olib.DEBUG_FULL]


@dataclasses.dataclass
class TransformParameters:
    """Parameters required for the calculations."""

    #: If telemetry cannot be determined, use existing information in the observation's header.
    allow_default: bool = False
    #: Detector aperture
    aperture: str | None = None
    #: FGS Boresite Adjustment Matrix (BAM) reference file path
    bam_ref: str | Path | None = None
    #: Parameters needed for CRDS queries
    crds_pars: dict | None = None
    #: Default quaternion to use if engineering is not available.
    default_quaternion: tuple | None = None
    #: The distortion reference which defines the V1 to detector transformation.
    distortion_ref: str | Path | None = None
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
    #: Target aperture
    target_aperture: str | None = None
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
    # Gather all parameters for CRDS usage.
    t_pars.crds_pars = model.get_crds_parameters()

    # Instrument details
    t_pars.aperture = model.meta.wcsinfo.aperture_name
    try:
        exp_type = model.meta.exposure.type.lower()
    except AttributeError:
        exp_type = None
    t_pars.exp_type = exp_type
    t_pars.target_aperture = model.meta.pointing.target_aperture

    # observation parameters
    t_pars.obsstart = model.meta.exposure.start_time
    t_pars.obsend = model.meta.exposure.end_time
    ephem = model.meta.ephemeris
    t_pars.velocity = (ephem.velocity_x, ephem.velocity_y, ephem.velocity_z)

    # Retrieve previously calculated orientation items only if they are currently not defined.
    if t_pars.default_quaternion is None:
        t_pars.default_quaternion = getattr(model.meta.pointing, "quaternion", None)
