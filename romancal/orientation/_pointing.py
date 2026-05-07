"""Retrieve orientation parameters from the engineering database"""

import logging
from collections import defaultdict, namedtuple
from copy import copy

import numpy as np
from astropy.table import Table
from astropy.time import Time, TimeDelta

from . import _lib as olib

from ..lib.engdb.engdb_lib import EngDB_Value
from ..lib.engdb.engdb_tools import engdb_service

__all__ = []

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LOGLEVELS = [logging.INFO, logging.DEBUG, olib.DEBUG_FULL]

# Mnemonics needed.
COARSE_MNEMONICS_QUATERNION_ECI = [f"SCF_AC_SDR_QBJ_{idx + 1}" for idx in range(4)]
COARSE_MNEMONICS_B2FGS_EST = [f"SCF_AC_EST_FGS_QBR_{idx + 1}" for idx in range(4)]
COARSE_MNEMONICS = COARSE_MNEMONICS_QUATERNION_ECI + COARSE_MNEMONICS_B2FGS_EST

# Pointing container
# Attributes are as follows. Except for the observation time, all values
# are retrieved from the engineering data.
#    fgs_q        : Quaternion representing orientation of the FGS frame relative to the Observatory frame
#    obstime      : Time the pointing information refers to.
#    q            : Quaternion of the FGS.
Pointing = namedtuple("Pointing", ["fgs_q", "obstime", "q"])
Pointing.__new__.__defaults__ = (None,) * 3

def all_pointings(mnemonics):
    """
    V1 of making pointings.

    Parameters
    ----------
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Returns
    -------
    pointings : [Pointing[,...]]
        List of pointings.
    """
    ordered = mnemonics_chronologically(mnemonics)
    pointings = mnemonics_to_pointings(ordered)

    return pointings


def first_pointing(mnemonics):
    """
    Return first pointing.

    Parameters
    ----------
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Returns
    -------
    pointing : Pointing
        First pointing.
    """
    pointings = all_pointings(mnemonics)
    return pointings[0]


def get_mnemonics(
    obsstart, obsend, tolerance, mnemonics_to_read=COARSE_MNEMONICS, service_kwargs=None
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

    mnemonics_to_read : (str[,...])
        The mnemonics to fetch.

    service_kwargs : dict or None
        Keyword arguments passed to `engdb_service` defining what
        engineering database service to use.

    Returns
    -------
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Raises
    ------
    ValueError
        Cannot retrieve engineering information.
    """
    if service_kwargs is None:
        service_kwargs = dict()
    try:
        engdb = engdb_service(**service_kwargs)
    except olib.EXPECTED_ERRORS as exception:
        raise ValueError(
            f"Cannot open engineering DB connection\nException: {exception}"
        ) from exception
    logger.info("Querying engineering DB: %s", engdb)

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
        except olib.EXPECTED_ERRORS as exception:
            logger.warning("Cannot retrieve %s from engineering.", mnemonic)
            logger.debug("Exception %s", exception)
            continue

        # If more than two points exist, throw off the bracket values.
        # Else, ensure the bracket values are within the allowed time.
        if len(mnemonics[mnemonic]) < 2:
            logger.warning(
                "Mnemonic %s has no telemetry within the observation time.", mnemonic
            )
            logger.warning(
                "Attempting to use bracket values within %s seconds", tolerance
            )

            mnemonics[mnemonic] = engdb.get_values(
                mnemonic,
                obsstart,
                obsend,
                time_format="mjd",
                include_obstime=True,
                include_bracket_values=True,
            )

            tolerance_mjd = TimeDelta(tolerance, format="sec")
            allowed_start = obsstart - tolerance_mjd
            allowed_end = obsend + tolerance_mjd
            allowed = [
                value
                for value in mnemonics[mnemonic]
                if allowed_start <= value.obstime <= allowed_end
            ]
            if not len(allowed):
                logger.warning(
                    "No telemetry exists for mnemonic {} within {} and {}".format(
                        mnemonic,
                        Time(allowed_start, format="mjd").isot,
                        Time(allowed_end, format="mjd").isot,
                    )
                )
            mnemonics[mnemonic] = allowed

    return mnemonics


def get_pointing(
    obsstart,
    obsend,
    mnemonics_to_read=COARSE_MNEMONICS,
    service_kwargs=None,
    tolerance=60,
    reduce_func=None,
):
    """
    Get telescope pointing engineering data.

    Parameters
    ----------
    obsstart, obsend : float
        MJD observation start/end times

    mnemonics_to_read : (str,[...])
        mnemonics to read.

    service_kwargs : dict or None
        Keyword arguments passed to `engdb_service` defining what
        engineering database service to use.

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
        service_kwargs=service_kwargs,
    )
    reduced = reduce_func(mnemonics)

    logger.log(olib.DEBUG_FULL, "Mnemonics found:")
    logger.log(olib.DEBUG_FULL, "%s", mnemonics)
    logger.info("Reduced set of pointings:")
    logger.info("%s", reduced)

    return reduced


def mnemonics_chronologically(mnemonics):
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

    Returns
    -------
    ordered : [(obstime, {mnemonic: value[,...]}[,...]]
        Time-ordered list of 2-tuple consisting of `(time, mnemonics)`.
    """
    # Collect all information by observation time and sort.
    by_obstime = defaultdict(dict)
    for mnemonic, values in mnemonics.items():
        if values is not None:
            for value in values:
                by_obstime[value.obstime][mnemonic] = value
    by_obstime = sorted(by_obstime.items())

    # Created the ordered matrix
    ordered = []
    last_obstime = {}
    for obstime, mnemonics_at_time in by_obstime:
        last_obstime.update(mnemonics_at_time)

        # Engineering data may be present, but all zeros.
        # Filter out this situation.
        values = [value.value for value in last_obstime.values()]
        if not any(values):
            continue

        ordered.append((obstime, copy(last_obstime)))

    return ordered


def mnemonics_chronologically_table(mnemonics):
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

    Returns
    -------
    ordered_by_time : `astropy.table.Table`
        Time-ordered mnemonic list with progressive values.
    """
    ordered = mnemonics_chronologically(mnemonics)

    names = list(mnemonics.keys())
    names = ["time", *names]
    time_idx = 0

    values = [[] for _ in names]

    for time, mnemonics_at_time in ordered:
        values[time_idx].append(time)
        for mnemonic in mnemonics_at_time:
            idx = names.index(mnemonic)
            values[idx].append(ordered[time][mnemonic].value)

    t = Table(values, names=names)

    return t


def mnemonics_to_pointings(ordered_mnemonics):
    """From a list of mnemonic values, create a list of Pointings

    Parameters
    ----------
    ordered_mnemonics : [(Time, {mnemonic: [value[,...]][,...]})[,...]]
        List of 2-tuples consisting of (Time, dict) where the dict are
        the mnemonics values at the associated time.

    Returns
    -------
    pointings : [Pointing[,...]]
        List of pointings
    """
    pointings = []
    for obstime, mnemonics_at_time in ordered_mnemonics:
        # Observatory orientation, required
        try:
            q = np.array(
                [mnemonics_at_time[m].value for m in COARSE_MNEMONICS_QUATERNION_ECI]
            )
        except KeyError as exception:
            raise ValueError(
                f"One or more quaternion mnemonics not in the telemetry {COARSE_MNEMONICS_QUATERNION_ECI}"
            ) from exception

        # B-frame to FGS-frame quaternion. Not required and very oddly has so many backups...
        fgs_q = None
        try:
            fgs_q = np.array(
                [mnemonics_at_time[m].value for m in COARSE_MNEMONICS_B2FGS_EST]
            )
        except KeyError:
            logger.warning(
                "One or more of the B-to-FGS quaternion mnemonics are not in the telementry %s",
                COARSE_MNEMONICS_B2FGS_EST,
            )

        pointing = Pointing(fgs_q=fgs_q, obstime=obstime, q=q)
        pointings.append(pointing)

    return pointings


def pointing_from_average(mnemonics):
    """
    Determine single pointing from average of available pointings.

    Parameters
    ----------
    mnemonics : {mnemonic: [value[,...]][,...]}
        The values for each pointing mnemonic.

    Returns
    -------
    pointing : Pointing
        Pointing from average.
    """
    # weed-out empty mnemonics
    valid_mnemonics = [
        mnemonic for mnemonic in mnemonics if mnemonics[mnemonic] is not None
    ]

    # Get average observation time.
    times = [
        eng_param.obstime.unix
        for mnemonic in valid_mnemonics
        for eng_param in mnemonics[mnemonic]
        if eng_param.obstime.unix != 0.0
    ]
    if len(times) > 0:
        obstime = Time(np.average(times), format="unix")
    else:
        raise ValueError("No valid mnemonics found.")

    # Get averages for all the mnemonics.
    mnemonic_averages = {}
    for mnemonic in valid_mnemonics:
        values = [eng_param.value for eng_param in mnemonics[mnemonic]]
        if np.allclose(values, 0.0):
            logger.warning(
                "Mnemonics %s is only zeros. Treating as undefined.", mnemonic
            )
        else:
            mnemonic_averages[mnemonic] = EngDB_Value(
                obstime=obstime, value=np.average(values)
            )

    pointing = mnemonics_to_pointings([(obstime, mnemonic_averages)])[0]

    # That's all folks
    return pointing
