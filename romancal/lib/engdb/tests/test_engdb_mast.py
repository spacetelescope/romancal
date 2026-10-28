"""Test the MAST Engineering interface"""

import logging
from pathlib import Path

import pytest

from astropy.table import Table
from astropy.time import Time
import numpy as np

from romancal.lib.engdb import engdb_mast
from romancal.lib.engdb.engdb_lib import EngDB_Value

# Configure logging
log = logging.getLogger(__name__)

# Test query
QUERY = ("scf_ac_sdr_qbj_1", "2026-10-01T00:00:00", "2026-10-01T00:05:00")


@pytest.mark.parametrize(
    "mnemonic, expected",
    [
        (None, "something"),
        ("scf_ac_sdr_qbj_1", 1),
        ("ope", "something"),
        ("junkfromspace", 0),
    ],
)
def test_get_meta(engdb, mnemonic, expected):
    """Test meta retrieval"""
    results = engdb.get_meta(search=mnemonic)
    n = results["Count"]
    assert n == len(results["TlmMnemonics"]), "No meta for mnemonics found"
    if expected == "something":
        assert n != 0, f"Unexpected database contents. Check state of database. Count: {n}, expected: {expected}"
    else:
        assert n == expected, f"Unexpected database contents. Check state of database. Count: {n}, expected: {expected}"


def test_get_records(engdb):
    """Test getting records"""
    records = engdb._get_records(*QUERY)
    assert len(records)
    assert np.allclose(records['EUValue'], 0.756946)
    assert (records['MJD'] >= (Time(QUERY[1]).mjd - 60.)).all()
    assert (records['MJD'] <= (Time(QUERY[2]).mjd + 60.)).all()


@pytest.mark.parametrize(
    "contents",
    [
        'TlmMnemonic',
        'EUValue',
        'ObsTime',
        'ReqSTime',
        'ReqETime',
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
    assert np.allclose(values, 0.7569, rtol=1.e-4)


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
    assert isinstance(item.value, float)


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
    assert isinstance(result.value[0], float)


@pytest.mark.parametrize(
    "envs, expected",
    [
        ({"MAST_AUTH_TOKEN": None, "MAST_API_TOKEN": None}, {}),
        (
            {"MAST_AUTH_TOKEN": "MAST_AUTH_TOKEN", "MAST_API_TOKEN": None},
            {"x-asb-auth": "MAST_AUTH_TOKEN"},
        ),
        (
            {"MAST_AUTH_TOKEN": None, "MAST_API_TOKEN": "MAST_API_TOKEN"},
            {"Authorization": "token MAST_API_TOKEN"},
        ),
        (
            {"MAST_AUTH_TOKEN": "MAST_AUTH_TOKEN", "MAST_API_TOKEN": "MAST_API_TOKEN"},
            {"x-asb-auth": "MAST_AUTH_TOKEN"},
        ),
    ],
)
def test_headers_env(envs, expected, monkeypatch):
    """Check that headers are being appropriately set from environment"""
    for env, value in envs.items():
        if value is None:
            monkeypatch.delenv(env, raising=False)
        else:
            monkeypatch.setenv(env, value)

    edb = engdb_mast.EngdbMast(check_aliveness=False)

    assert edb._metareq.headers == expected


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        ({}, {}),
        ({"token": "token"}, {"Authorization": "token token"}),
        ({"token": "token", "rsdp_auth": True}, {"x-asb-auth": "token"}),
    ],
)
def test_headers_kwargs(kwargs, expected, monkeypatch):
    """Check that headers are being appropriately set"""
    # Ensure tokens do not get confused by external settings
    monkeypatch.delenv("MAST_AUTH_TOKEN", raising=False)
    monkeypatch.delenv("MAST_API_TOKEN", raising=False)

    edb = engdb_mast.EngdbMast(check_aliveness=False, **kwargs)

    assert edb._metareq.headers == expected


def test_negative_aliveness():
    """Ensure failure occurs with a bad url"""
    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast(
            eng_base_url="https://127.0.0.1/_engdb_mast_test",
            token="dummytoken",  # noqa: S106
        )


@pytest.mark.parametrize(
    "env, var, value, default",
    [
        ("ENG_BASE_URL", "eng_base_url", None, engdb_mast.MAST_BASE_URL + "/"),
        ("ENG_BASE_URL", "eng_base_url", "something/", None),
        ("ENG_DATA_URI", "data_uri", None, engdb_mast.DATA_URI),
        ("ENG_DATA_URI", "data_uri", "something", None),
        ("ENG_META_URI", "meta_uri", None, engdb_mast.META_URI),
        ("ENG_META_URI", "meta_uri", "something", None),
    ],
)
def test_uri_from_environment(env, var, value, default, monkeypatch):
    """Test settings from environment variables"""
    if value is None:
        monkeypatch.delenv(env, raising=False)
        value = default
    else:
        monkeypatch.setenv(env, value)

    s = engdb_mast.EngdbMast(check_aliveness=False)

    assert getattr(s, var) == value


# ######################
# Fixtures and utilities
# ######################
@pytest.fixture(scope="module")
def engdb():
    """Open a connection"""
    try:
        engdb = engdb_mast.EngdbMast()
    except RuntimeError as exception:
        pytest.skip(f"Live MAST Engineering Service not available: {exception}")
    return engdb
