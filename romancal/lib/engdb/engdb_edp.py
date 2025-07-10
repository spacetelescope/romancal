"""Access the Roman Engineering Mnemonic Database through the direct EDP interface."""

import json
import logging
from ast import literal_eval
from os import getenv
from pathlib import Path
from shutil import copy2

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time
from requests.adapters import HTTPAdapter, Retry

from edp.mnemonics_reader import MnemonicsReader

from .engdb_lib import (
    FORCE_STATUSES,
    RETRIES,
    TIMEOUT,
    EngDB_Value,
    EngdbABC,
    mnemonic_data_fname,
)

__all__ = ["EngdbEDP"]

# Default MAST info.
MAST_BASE_URL = "https://masttest.stsci.edu"
API_URI = "edp/api/v0.1/mnemonics/fqa/roman/data"

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class EngdbEDP(EngdbABC):
    """
    Access the Roman Engineering Database through the direct EDP interface.

    The EDP interface is maintained by STScI MAST.

    Two items are required. An `edp.keytab` for authentication to the EDP servers.

    Environment variable MSSQL_DRIVER needs to point to the Microsoft SQL driver library.

    Parameters
    ----------
    environment : ['dev', 'test', 'int', 'ops']
        The operating environment currently working under.

    path_to_cc : Path-like
        Path to the kerberos keytab required to access the database.

    **service_kwargs : dict
        Service-specific keyword arguments that are not relevant to this implementation
        of EngdbABC.

    Raises
    ------
    RuntimeError
        Any and all failures with connecting with the MAST server.
    """
    #: The end time of the last query.
    endtime = None

    #: The MnemonicsReader instance
    mr = None

    #: Operating environment
    ops_env = None

    #: Path to the keytab
    path_to_cc: None

    #: The results of the last query.
    response = None

    #: Number of retries to attempt to contact the service
    retries = RETRIES

    #: The start time of the last query.
    starttime = None

    #: Network timeout when communicating with the service
    timeout = TIMEOUT

    def __init__(self, environment, path_to_cc, **service_kwargs):
        logger.debug("kwargs not used by this service: %s", service_kwargs)

        try:
            self.configure(environment, path_to_cc)
        except FileNotFoundError as exception:
            raise RuntimeError(f'Cannot instantiate EDP service with environment: {environment}, path_to_cc: {path_to_cc}') from exception

        # Check for basic aliveness.
        self.mr.get_mnemonic_metadata()

    def configure(self, ops_env, path_to_cc):
        """
        Configure from parameters and environment.

        Parameters
        ----------
        environment : ['dev', 'test', 'int', 'ops']
            The operating environment currently working under.

        path_to_cc : Path-like
            Path to the kerberos keytab required to access the database.
        """
        self.ops_env = ops_env
        self.path_to_cc = path_to_cc

        # Create the mnemonic reader
        self.mr = MnemonicsReader(mission='roman', env=ops_env, access_type='FqA', path_to_cc=path_to_cc)

        # Get various timeout parameters
        self.retries = getenv("ENG_RETRIES", RETRIES)
        self.timeout = getenv("ENG_TIMEOUT", TIMEOUT)

    def get_meta(self, *kwargs):
        """
        Get the mnemonics meta info.

        The MAST interface does not provide any meta.
        """
        raise NotImplementedError(
            "MAST Engineering AUI does not provide a meta service"
        )

    def get_values(
        self,
        mnemonic,
        starttime,
        endtime,
        time_format=None,
        include_obstime=False,
        include_bracket_values=False,
        zip_results=True,
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic : str
            The engineering mnemonic to retrieve

        starttime : str or `astropy.time.Time`
            The, inclusive, start time to retrieve from.

        endtime : str or `astropy.time.Time`
            The, inclusive, end time to retrieve from.

        time_format : str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        include_obstime : bool
            If `True`, the return values will include observation
            time as `astropy.time.Time`. See `zip_results` for further details.

        include_bracket_values : bool
            The DB service, by default, returns the bracketing
            values outside of the requested time. If `True`, include
            these values.

        zip_results : bool
            If `True` and `include_obstime` is `True`, the return values
            will be a list of 2-tuples. If false, the return will
            be a single 2-tuple, where each element is a list.

        Returns
        -------
        values : [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
            Returns the list of values. See `include_obstime` and `zip` for modifications.
        """
        if not isinstance(starttime, Time):
            starttime = Time(starttime, format=time_format)
        if not isinstance(endtime, Time):
            endtime = Time(endtime, format=time_format)

        records = self._get_records(
            mnemonic=mnemonic,
            starttime=starttime,
            endtime=endtime,
            time_format=time_format,
        )
        if len(records) == 0:
            return list()

        # If desired, remove bracket or outside of timeframe entries.
        if not include_bracket_values:
            selection = np.logical_and(
                records["MJD"] >= starttime.mjd, records["MJD"] <= endtime.mjd
            )
            records = records[selection]

        # Reformat to the desired list formatting.
        results = _ValueCollection(
            include_obstime=include_obstime, zip_results=zip_results
        )
        values = records["EUValue"]
        obstimes = Time(records["MJD"], format="mjd")
        for obstime, value in zip(obstimes, values, strict=False):
            results.append(obstime, value)

        return results.collection

    def set_session(self):
        """Set up HTTP session."""
        self._req = requests.Request(
            method="GET",
            url=self.base_url + API_URI,
            headers={"Authorization": f"token {self.token}"},
        )

        s = requests.Session()
        retries = Retry(
            total=self.retries,
            backoff_factor=1.0,
            status_forcelist=FORCE_STATUSES,
            raise_on_status=True,
        )
        s.mount("https://", HTTPAdapter(max_retries=retries))

        self._session = s

    def _get_records(
        self, mnemonic, starttime, endtime, time_format=None, **other_kwargs
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic : str
            The engineering mnemonic to retrieve

        starttime : str or astropy.time.Time
            The, inclusive, start time to retrieve from.

        endtime : str or astropy.time.Time
            The, inclusive, end time to retrieve from.

        time_format : str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        **other_kwargs : dict
            Keyword arguments not relevant to this implementation.

        Returns
        -------
        records : `astropy.Table`
            Returns the resulting table.

        Notes
        -----
        The engineering service always returns the bracketing entries
        before and after the requested time range.
        """
        if not isinstance(starttime, Time):
            starttime = Time(starttime, format=time_format)
        if not isinstance(endtime, Time):
            endtime = Time(endtime, format=time_format)
        self.starttime = starttime
        self.endtime = endtime

        # Make the query.
        response = self.mr.get_mnemonic_data(mnemonic, starttime.isot, endtime.isot)
        self.response = response.copy()

        # Convert to table.
        data = response["Data"]
        del response["Data"]
        table = Table(rows=data, meta=response)

        # Create a column MJD that has the MJD version of the data
        if len(table):
            obstime = Time(table["ObsTime"])
            table["MJD"] = obstime.mjd

        return table

    def __repr__(self):
        """What am I"""
        mr = self.mr
        db_config = mr.db_config.copy()
        del db_config['ad_name']
        repr = f"{self.__class__.__name__}(config={db_config}, mission='{mr.mission}', env='{mr.env}')"
        return repr


class _ValueCollection:
    """
    Engineering Value Collection.

    Parameters
    ----------
    include_obstime : bool
        If `True`, the return values will include observation
        time as `astropy.time.Time`. See `zip_results` for further details.

    zip_results : bool
        If `True` and `include_obstime` is `True`, the return values
        will be a list of 2-tuples. If false, the return will
        be a single 2-tuple, where each element is a list.

    Attributes
    ----------
    collection : [value, ...] or [(obstime, value), ...] or ([obstime,...], [value, ...])
        Returns the list of values.
        See `include_obstime` and `zip_results` for modifications.
    """

    def __init__(self, include_obstime=False, zip_results=True):
        self._include_obstime = include_obstime
        self._zip_results = zip_results
        if zip_results:
            self.collection = []
        else:
            self.collection = EngDB_Value([], [])

    def append(self, obstime, value):
        """
        Append value to collection.

        Parameters
        ----------
        obstime : `astropy.time.Time`
            Observation time as returned from the engineering.

        value : numeric
            Value from DB.
        """
        # Make all the times readable
        obstime.format = "isot"

        # Append
        if self._include_obstime:
            if self._zip_results:
                self.collection.append(EngDB_Value(obstime, value))
            else:
                self.collection.obstime.append(obstime)
                self.collection.value.append(value)
        else:
            self.collection.append(value)
