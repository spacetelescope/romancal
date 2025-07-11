"""V1 Calculation based on time and engineering database info."""

import logging
from collections import defaultdict

import roman_datamodels as rdm
from astropy.table import Table

from . import set_telescope_pointing as stp

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ["v1_calculate_from_models", "v1_calculate_over_time"]


def v1_calculate_from_models(sources, **calc_wcs_from_time_kwargs):
    """
    Calculate V1 over the time period for the given models.

    Returns a table of V1 pointings for all input models.
    The table has the following columns:

    * source (`roman_datamodels.datamodels.RomanDataModel`): The model.
    * obstime (`astropy.time.Time`): The observation time.
    * v1 (float, float, float): 3-tuple or RA, dec, and position angle.

    Parameters
    ----------
    sources : [File-like or roman_datamodels.datamodels.RomanDataModel[...]]
        The datamodels to get timings other header parameters from.

    **calc_wcs_from_time_kwargs : dict
        Keyword arguments to pass to ``calc_wcs_from_time``.

    Returns
    -------
    v1_table : `astropy.table.Table`
        Table of V1 pointing.
    """
    # Initialize structures.
    v1_dict = defaultdict(list)

    # Calculate V1 for all sources.
    for source in sources:
        with rdm.open(source) as model:
            t_pars = stp.t_pars_from_model(model, **calc_wcs_from_time_kwargs)
            obstimes, _, vinfos = stp.calc_wcs_over_time(
                t_pars.obsstart, t_pars.obsend, t_pars
            )

        sources = [source] * len(obstimes)
        v1_dict["source"] += sources
        v1_dict["obstime"] += obstimes
        v1_dict["v1"] += vinfos

    # Format and return.
    v1_table = Table(v1_dict, meta=t_pars.as_reprdict())
    return v1_table


def v1_calculate_over_time(obsstart, obsend, **calc_wcs_from_time_kwargs):
    """
    Calculate V1 over the given time period.

    Returns a table of all V1 pointings that can be retrieved from the engineering database
    that exist between, inclusively, the start and end times.

    The table has the following columns:

    * source (str): The string "time range".
    * obstime (`astropy.time.Time`): The observation time.
    * v1 (float, float, float): 3-tuple or RA, dec, and position angle.

    Parameters
    ----------
    obsstart, obsend : float
        The MJD start and end time to search for pointings.

    **calc_wcs_from_time_kwargs : dict
        Keyword arguments to pass to ``calc_wcs_from_time``.

    Returns
    -------
    v1_table : `astropy.table.Table`
        Table of V1 pointing.
    """
    # Initialize structures.
    t_pars = stp.TransformParameters(**calc_wcs_from_time_kwargs)
    if len(t_pars.aperture) == 0:
        t_pars.aperture = "WFI_CEN"

    # Calculate V1 for all sources.
    obstimes, _, vinfos = stp.calc_wcs_over_time(obsstart, obsend, t_pars)
    v1_dict = {}
    v1_dict["source"] = ["time range"] * len(obstimes)
    v1_dict["obstime"] = obstimes
    v1_dict["v1"] = vinfos

    # Format and return.
    v1_table = Table(v1_dict, meta=t_pars.as_reprdict())
    return v1_table


def simplify_table(v1_table):
    """
    Convert pure object-based table to ASCII/Human-friendly.

    The tables as produced by the `v1_calculate` functions use native objects.
    For instance, the "obstime" column contains `astropy.time.Time` objects and
    "v1" is the `romancal.orientation.set_telescope_pointing.WCSREF` object.

    This routine converts such objects to strings or Python-native builtin objects.

    Parameters
    ----------
    v1_table : `astropy.table.Table`
        V1 table as produced by ``v1_calculate`` functions.

    Returns
    -------
    formatted : `astropy.table.Table`
        Reformatted table.
    """
    source_formatted = [str(v) for v in v1_table["source"]]
    obstime_formatted = v1_table["obstime"].isot
    ras, decs, pa_v3s = list(map(list, zip(*v1_table["v1"], strict=False)))

    formatted = Table(
        [source_formatted, obstime_formatted, ras, decs, pa_v3s],
        names=("source", "obstime", "ra", "dec", "pa_v3"),
        meta=v1_table.meta,
    )
    return formatted
