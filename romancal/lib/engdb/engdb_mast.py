"""Access the Roman Engineering Mnemonic Database through MAST."""

import logging
from ast import literal_eval
from os import getenv

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time
from requests.adapters import HTTPAdapter, Retry

from .engdb_lib import (
    FORCE_STATUSES,
    RETRIES,
    TIMEOUT,
    EngdbABC,
    ValueCollection,
)

__all__ = ["EngdbMast"]

# Default MAST info.
MAST_BASE_URL = "https://mast.stsci.edu"
DATA_URI = "edp/api/v0.1/mnemonics/spa/roman/data"
META_URI = "edp/api/v0.1/mnemonics/spa/roman/metadata"

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class EngdbMast(EngdbABC):
    """
    Access the Roman Engineering Database through MAST.

    Parameters
    ----------
    base_url : str
        The base url for the engineering RESTful service. If not defined,
        the environmental variable ENG_BASE_URL is queried. Otherwise
        the default MAST website is used.

    token : str or None
        The MAST access token. If not defined, the environmental variable
        MAST_API_TOKEN is queried. A token is required.
        For more information, see https://auth.mast.stsci.edu/

    **service_kwargs : dict
        Service-specific keyword arguments that are not relevant to this implementation
        of EngdbABC.

    Raises
    ------
    RuntimeError
        Any and all failures with connecting with the MAST server.
    """

    #: The base URL for the engineering service.
    base_url = None

    #: The end time of the last query.
    endtime = None

    #: The results of the last query.
    response = None

    #: Number of retries to attempt to contact the service
    retries = RETRIES

    #: The start time of the last query.
    starttime = None

    #: Network timeout when communicating with the service
    timeout = TIMEOUT

    #: MAST Token
    token = None

    def __init__(self, base_url=None, token=None, **service_kwargs):
        logger.debug("kwargs not used by this service: %s", service_kwargs)

        self.configure(base_url=base_url, token=token)

        # Check for basic aliveness.
        try:
            resp = requests.get(self.base_url + "edp/", timeout=self.timeout)
        except requests.exceptions.ConnectionError as exception:
            raise RuntimeError(
                f"MAST url: {self.base_url} is unreachable."
            ) from exception
        if resp.status_code != 200:
            raise RuntimeError(
                f"MAST url: {self.base_url} is not available. "
                f"Returned HTTPS status {resp.status_code}"
            )

        # Basics are covered. Finalize initialization.
        self.set_session()

    def configure(self, base_url=None, token=None):
        """
        Configure from parameters and environment.

        Parameters
        ----------
        base_url : str
            The base url for the engineering RESTful service. If not defined,
            the environmental variable ENG_BASE_URL is queried. Otherwise
            the default MAST website is used.

        token : str or None
            The MAST access token. If not defined, the environmental variable
            MAST_API_TOKEN is queried. A token is required.
            For more information, see 'https://auth.mast.stsci.edu/'
        """
        # Determine the database to use
        if base_url is None:
            base_url = getenv("ENG_BASE_URL", MAST_BASE_URL)
        if base_url[-1] != "/":
            base_url += "/"
        self.base_url = base_url

        # Get the token
        if token is None:
            token = getenv("MAST_API_TOKEN", None)
        self.token = token

        # Get various timeout parameters
        self.retries = int(getenv("ENG_RETRIES", RETRIES))
        self.timeout = int(getenv("ENG_TIMEOUT", TIMEOUT))

    def get_meta(self, search=None):
        """
        Get the mnemonics meta info.

        Parameters
        ----------
        search : str or None
            A partial, or full, mnemonic specification.
            If None, meta for all available mnemonics are returned

        Returns
        -------
        meta : ???
            The meta information
        """

        # Make the request
        if search is not None:
            self._metareq.params = {
                "mnemonic": search,
            }
        prepped = self._session.prepare_request(self._metareq)
        settings = self._session.merge_environment_settings(
            prepped.url, {}, None, None, None
        )
        logger.debug("Query: %s", prepped.url)
        self.metaresponse = self._session.send(prepped, timeout=self.timeout, **settings)
        self.metaresponse.raise_for_status()
        logger.debug("Response: %s", self.metaresponse)
        logger.debug("Response test: %s", self.metaresponse.text)

        # Leave as dictionary.
        results = literal_eval(self.metaresponse.text)
        return results


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
        results = ValueCollection(
            include_obstime=include_obstime, zip_results=zip_results
        )
        values = records["EUValue"]
        obstimes = Time(records["MJD"], format="mjd")
        for obstime, value in zip(obstimes, values, strict=False):
            results.append(obstime, value)

        return results.collection

    def set_session(self):
        """Set up HTTP session."""
        self._datareq = requests.Request(
            method="GET",
            url=self.base_url + DATA_URI,
            headers={"Authorization": f"token {self.token}"},
        )

        self._metareq = requests.Request(
            method="GET",
            url=self.base_url + META_URI,
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

        # Make the request
        mnemonic = mnemonic.strip()
        mnemonic = mnemonic.upper()
        starttime_fmt = starttime.strftime("%Y-%m-%dT%H:%M:%S")
        endtime_fmt = endtime.strftime("%Y-%m-%dT%H:%M:%S")
        self._datareq.params = {
            "mnemonic": mnemonic,
            "s_time": starttime_fmt,
            "e_time": endtime_fmt,
        }
        prepped = self._session.prepare_request(self._datareq)
        settings = self._session.merge_environment_settings(
            prepped.url, {}, None, None, None
        )
        logger.debug("Query: %s", prepped.url)
        self.response = self._session.send(prepped, timeout=self.timeout, **settings)
        self.response.raise_for_status()
        logger.debug("Response: %s", self.response)
        logger.debug("Response test: %s", self.response.text)

        # Convert to table.
        response = literal_eval(self.response.text)
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
        repr = f"{self.__class__.__name__}(base_url='{self.base_url}')"
        return repr
