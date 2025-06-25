"""Test the MAST Engineering interface"""
from pathlib import Path
import pytest
import requests

from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values

from romancal.lib.engdb_lib import EngDB_Value
from romancal.lib import engdb_mast

# Test query
QUERY = ('ope_scf_dir', '2027-02-23T01:00:00', '2027-02-23T01:00:05')

# Expected return from query
EXPECTED_RESPONSE = ('{"TlmMnemonic":"OPE_SCF_DIR",'
                     '"ReqSTime":"2027-02-23T01:00:00.000",'
                     '"ReqETime":"2027-02-23T01:00:05.000",'
                     '"Count":7,"AllPoints":1,"Data":['
                     '{"ObsTime":"2027-02-23T00:59:59.054","EUValue":"SCFA"},'
                     '{"ObsTime":"2027-02-23T01:00:00.052","EUValue":"SCFA"},'
                     '{"ObsTime":"2027-02-23T01:00:01.052","EUValue":"SCFA"},'
                     '{"ObsTime":"2027-02-23T01:00:02.162","EUValue":"SCFA"},'
                     '{"ObsTime":"2027-02-23T01:00:03.161","EUValue":"SCFA"},'
                     '{"ObsTime":"2027-02-23T01:00:04.161","EUValue":"SCFA"},'
                     '{"ObsTime":"2027-02-23T01:00:05.163","EUValue":"SCFA"}]}')
EXPECTED_RECORDS = Table.read(Path(__file__).parent / 'data' / 'test_records_expected.ecsv', format='ascii.ecsv')


@pytest.fixture(scope='module')
def is_alive():
    """Check if the MAST portal is accessible
    """
    is_alive = False
    try:
        r = requests.get(engdb_mast.MAST_BASE_URL)
        is_alive = (r.status_code == requests.codes.ok)
    except Exception:
        pass
    if not is_alive:
        pytest.skip(f'MAST url {engdb_mast.MAST_BASE_URL} not available. Skipping.')


@pytest.fixture(scope='module')
def engdb():
    """Open a connection"""
    try:
        engdb = engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f'Live MAST Engineering Service not available: {exception}')
    return engdb


def test_aliveness(is_alive):
    """Check connection creation

    Failure is any failure from instantiation.
    """
    engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL, token='dummytoken')


def test_get_records(engdb):
    """Test getting records"""
    records = engdb._get_records(*QUERY)
    assert engdb.response.text == EXPECTED_RESPONSE
    assert report_diff_values(records, EXPECTED_RECORDS)


@pytest.mark.parametrize(
    'pars, expected',
    [
        ({}, ['SCFA'] * 5 ),
        ({'include_obstime': True},
         [EngDB_Value(obstime=Time(61459.04166726852, format='mjd'), value='SCFA'),
          EngDB_Value(obstime=Time(61459.041678842594, format='mjd'), value='SCFA'),
          EngDB_Value(obstime=Time(61459.04169168982, format='mjd'), value='SCFA'),
          EngDB_Value(obstime=Time(61459.04170325232, format='mjd'), value='SCFA'),
          EngDB_Value(obstime=Time(61459.04171482639, format='mjd'), value='SCFA')]),
        ({'include_obstime': True, 'zip_results': False}, EngDB_Value(
            obstime=[
                Time(61459.04166726852, format='mjd'),
                Time(61459.041678842594, format='mjd'),
                Time(61459.04169168982, format='mjd'),
                Time(61459.04170325232, format='mjd'),
                Time(61459.04171482639, format='mjd')],
            value=['SCFA'] * 5)),
        ({'include_bracket_values': True},
         ['SCFA'] * 7)
    ])
def test_get_values(engdb, pars, expected):
    values = engdb.get_values(*QUERY, **pars)
    assert values == expected


def test_negative_aliveness():
    """Ensure failure occurs with a bad url"""
    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast(base_url='https://127.0.0.1/_engdb_mast_test', token='dummytoken')
