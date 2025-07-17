"""Test the MAST Engineering interface"""

import logging
import os
from pathlib import Path

import pytest
from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values

from romancal.lib.engdb import engdb_edp
from romancal.lib.engdb.engdb_lib import EngDB_Value

pytest.importorskip("edp.mnemonics_reader")

# Configure logging
log = logging.getLogger(__name__)

# Test query
QUERY = ("ope_scf_dir", "2027-02-23T01:00:00", "2027-02-23T01:00:05")

# Expected return from query
EXPECTED_RESPONSE = {
    "TlmMnemonic": "ope_scf_dir",
    "ReqSTime": "2027-02-23T01:00:00.000",
    "ReqETime": "2027-02-23T01:00:05.000",
    "Count": 7,
    "AllPoints": 1,
    "Data": [
        {"ObsTime": "2027-02-23T00:59:59.054", "EUValue": "SCFA"},
        {"ObsTime": "2027-02-23T01:00:00.052", "EUValue": "SCFA"},
        {"ObsTime": "2027-02-23T01:00:01.052", "EUValue": "SCFA"},
        {"ObsTime": "2027-02-23T01:00:02.162", "EUValue": "SCFA"},
        {"ObsTime": "2027-02-23T01:00:03.161", "EUValue": "SCFA"},
        {"ObsTime": "2027-02-23T01:00:04.161", "EUValue": "SCFA"},
        {"ObsTime": "2027-02-23T01:00:05.163", "EUValue": "SCFA"},
    ],
}
EXPECTED_RECORDS = Table.read(
    Path(__file__).parent / "data" / "test_records_expected.ecsv", format="ascii.ecsv"
)


def test_get_records(engdb):
    """Test getting records"""
    records = engdb._get_records(*QUERY)
    assert engdb.response == EXPECTED_RESPONSE
    assert report_diff_values(records, EXPECTED_RECORDS)


@pytest.mark.parametrize(
    "pars, expected",
    [
        ({}, ["SCFA"] * 5),
        (
            {"include_obstime": True},
            [
                EngDB_Value(
                    obstime=Time(61459.04166726852, format="mjd"), value="SCFA"
                ),
                EngDB_Value(
                    obstime=Time(61459.041678842594, format="mjd"), value="SCFA"
                ),
                EngDB_Value(
                    obstime=Time(61459.04169168982, format="mjd"), value="SCFA"
                ),
                EngDB_Value(
                    obstime=Time(61459.04170325232, format="mjd"), value="SCFA"
                ),
                EngDB_Value(
                    obstime=Time(61459.04171482639, format="mjd"), value="SCFA"
                ),
            ],
        ),
        (
            {"include_obstime": True, "zip_results": False},
            EngDB_Value(
                obstime=[
                    Time(61459.04166726852, format="mjd"),
                    Time(61459.041678842594, format="mjd"),
                    Time(61459.04169168982, format="mjd"),
                    Time(61459.04170325232, format="mjd"),
                    Time(61459.04171482639, format="mjd"),
                ],
                value=["SCFA"] * 5,
            ),
        ),
        ({"include_bracket_values": True}, ["SCFA"] * 7),
    ],
)
def test_get_values(engdb, pars, expected):
    values = engdb.get_values(*QUERY, **pars)
    assert values == expected


@pytest.fixture(scope="module")
def engdb():
    """Open a connection"""
    try:
        engdb = engdb_edp.EngdbEDP("test", os.environ.get("KRB5CCPATH", None))
    except RuntimeError as exception:
        pytest.skip(f"Live Engineering Service not available: {exception}")
    return engdb
