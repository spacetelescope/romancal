"""Test the MAST Engineering interface"""

import logging
from pathlib import Path

import pytest
import requests
from astropy.table import Table
from astropy.time import Time

from romancal.lib.engdb import engdb_mast
from romancal.lib.engdb.engdb_lib import EngDB_Value

# Configure logging
log = logging.getLogger(__name__)

# Test query
QUERY = ("ope_scf_dir", "2027-02-1T00:00:00", "2027-02-28T23:00:00")

# Expected return from query
EXPECTED_RESPONSE = (
    '{"TlmMnemonic":"OPE_SCF_DIR",'
    '"ReqSTime":"2027-02-23T01:00:00.000",'
    '"ReqETime":"2027-02-23T01:00:05.000",'
    '"Count":7,"AllPoints":1,"Data":['
    '{"ObsTime":"2027-02-23T00:59:59.054","EUValue":"SCFA"},'
    '{"ObsTime":"2027-02-23T01:00:00.052","EUValue":"SCFA"},'
    '{"ObsTime":"2027-02-23T01:00:01.052","EUValue":"SCFA"},'
    '{"ObsTime":"2027-02-23T01:00:02.162","EUValue":"SCFA"},'
    '{"ObsTime":"2027-02-23T01:00:03.161","EUValue":"SCFA"},'
    '{"ObsTime":"2027-02-23T01:00:04.161","EUValue":"SCFA"},'
    '{"ObsTime":"2027-02-23T01:00:05.163","EUValue":"SCFA"}]}'
)
EXPECTED_RECORDS = Table.read(
    Path(__file__).parent / "data" / "test_records_expected.ecsv", format="ascii.ecsv"
)


@pytest.mark.parametrize('mnemonic, expected', [
    (None, 'something'),
    ('ope_scf_dir', 1),
    ('ope', 'something'),
    ('junkfromspace', 0),
])
def test_get_meta(engdb, mnemonic, expected):
    """Test meta retrieval"""
    results = engdb.get_meta(search=mnemonic)
    n = results['Count']
    assert n == len(results['TlmMnemonics'])
    if expected == 'something':
        if n == 0:
            pytest.xfail(reason=f'Unexpected database contents. Check state of database. Count: {n}, expected: {expected}')
    elif n != expected:
        pytest.xfail(reason=f'Unexpected database contents. Check state of database. Count: {n}, expected: {expected}')


def test_get_records(engdb):
    """Test getting records"""
    records = engdb._get_records(*QUERY)
    assert "SCFA" in records["EUValue"]


@pytest.mark.parametrize(
    "contents",
    [
        '"TlmMnemonic":"OPE_SCF_DIR"',
        '"EUValue":"SCFA"',
    ],
)
def test_get_records_response(engdb, contents):
    """Test getting records"""
    _ = engdb._get_records(*QUERY)
    assert contents in engdb.response.text

def test_get_value_justvalues(engdb):
    """Test just getting values"""
    values = engdb.get_values(*QUERY, include_bracket_values=True)
    assert len(values) > 1
    assert "SCFA" in values


def test_get_values_obstimes(engdb):
    """Get values with the observation times"""
    result = engdb.get_values(
        *QUERY, include_bracket_values=True, include_obstime=True, zip_results=True
    )
    assert isinstance(result, list)
    assert len(result) > 1
    item = result[0]
    assert isinstance(item, EngDB_Value)
    assert isinstance(item.obstime, Time)
    assert isinstance(item.value, str)


def test_get_values_nozip(engdb):
    """Get values with the observation times"""
    result = engdb.get_values(
        *QUERY, include_bracket_values=True, include_obstime=True, zip_results=False
    )
    assert isinstance(result, EngDB_Value)
    assert len(result) > 1
    assert isinstance(result.obstime, list)
    assert isinstance(result.value, list)
    assert isinstance(result.obstime[0], Time)
    assert isinstance(result.value[0], str)


def test_negative_aliveness():
    """Ensure failure occurs with a bad url"""
    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast(
            base_url="https://127.0.0.1/_engdb_mast_test",
            token="dummytoken",  # noqa: S106
        )


# ######################
# Fixtures and utilities
# ######################
@pytest.fixture(scope="module")
def engdb():
    """Open a connection"""
    try:
        engdb = engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f"Live MAST Engineering Service not available: {exception}")
    return engdb
