"""
Test suite for engdb_tools

Notes
-----
This file has been named specifically so it is not
automatically found by py.test. This is because, to test,
a connection to the internal engineering service is needed,
which is generally not available.
"""

import logging
import os

import pytest
from astropy.time import Time

from romancal.lib.engdb import engdb_tools

# Configure logging
log = logging.getLogger(__name__)

GOOD_MNEMONIC = "OPE_SCF_DIR"
GOOD_STARTTIME = "2027-02-23T01:00:00"
GOOD_ENDTIME = "2027-02-23T01:01:00"

SHORT_STARTTIME = "2027-02-23T01:00:30"

BAD_MNEMONIC = "No_Such_MNEMONIC"
NODATA_STARTIME = "2014-01-01"
NODATA_ENDTIME = "2014-01-02"


def test_environmental_bad(environ_jail):
    alternate = "https://google.com/"
    did_except = False
    os.environ["ENG_BASE_URL"] = alternate
    try:
        engdb = engdb_tools.engdb_service('mast')
    except Exception:
        did_except = True
    assert did_except, f"DB connection falsely created for {engdb.base_url}"


def test_basic(engdb):
    assert engdb._get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)


def test_bad_server():
    with pytest.raises(RuntimeError):
        engdb_tools.engdb_service('junk_service')


def test_values(engdb):
    records = engdb._get_records(GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME)
    assert len(records) == 2
    values = engdb.get_values(GOOD_MNEMONIC, GOOD_STARTTIME, SHORT_STARTTIME)
    assert len(values) == 29
    assert values[0] == "SCFA"


def test_values_with_bracket(engdb):
    records = engdb._get_records(GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME)
    assert len(records) == 2
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME, include_bracket_values=True
    )
    assert len(values) == 2
    assert values[1] == "SCFA"


def test_values_with_time(engdb):
    values = engdb.get_values(
        GOOD_MNEMONIC, GOOD_STARTTIME, SHORT_STARTTIME, include_obstime=True
    )
    assert len(values) >= 1
    assert isinstance(values[0], tuple)
    assert isinstance(values[0].obstime, Time)


def test_novalues(engdb):
    values = engdb.get_values(GOOD_MNEMONIC, NODATA_STARTIME, NODATA_ENDTIME)
    assert len(values) == 0


def test_unzip(engdb):
    """Test forunzipped versions of content"""
    values = engdb.get_values(
        GOOD_MNEMONIC,
        SHORT_STARTTIME,
        SHORT_STARTTIME,
        include_obstime=True,
        zip_results=False,
    )
    assert isinstance(values, tuple)
    assert len(values.obstime) == len(values.value)


# ########
# Fixtures
# ########
@pytest.fixture(scope='module')
def engdb():
    """Setup a service"""
    try:
        engdb = engdb_tools.engdb_service()
    except RuntimeError as exception:
        pytest.skip(f"Engineering database unvailable: {exception}.")
    yield engdb
