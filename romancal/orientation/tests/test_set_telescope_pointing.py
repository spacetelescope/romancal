"""
Test suite for set_telescope_pointing
"""

import dataclasses
import logging
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pysiaf")

import asdf
import roman_datamodels as rdm
from astropy.time import Time

from romancal.lib import (
    engdb_mast,
    engdb_tools
)
from romancal.orientation import set_telescope_pointing as stp

# Ensure that `set_telescope_pointing` logs.
stp.logger.setLevel(logging.DEBUG)
stp.logger.addHandler(logging.StreamHandler())

# Setup mock engineering service
STARTTIME = Time("2027-03-23T19:20:40", format="isot")
ENDTIME = Time("2027-03-23T19:21:36", format="isot")
BADSTARTTIME = Time("2020-02-02T02:02:02", format="isot")
BADENDTIME = Time("2020-02-02T02:12:02", format="isot")

# Header defaults
TARG_RA = 270.0
TARG_DEC = 66.01

# Get the mock databases
DATA_PATH = Path(__file__).parent / "data"

Q_EXPECTED = np.array([-0.69018802, 0.12195182, -0.695103, 0.15999998])
OBSTIME_EXPECTED = Time(1805829668.0276294, format="unix")

# Meta attributes for test comparisons
METAS_EQUALITY = [
    "meta.exposure.engineering_quality",
]
METAS_ISCLOSE = [
    "meta.wcsinfo.roll_ref",
    "meta.wcsinfo.v2_ref",
    "meta.wcsinfo.v3_ref",
    "meta.wcsinfo.ra_ref",
    "meta.wcsinfo.dec_ref",
    "meta.pointing.ra_v1",
    "meta.pointing.dec_v1",
    "meta.pointing.pa_v3",
]


def test_add_wcs_default(science_raw_model, tmp_path):
    """Handle when no pointing exists and the default is used."""
    m = science_raw_model
    m.meta.exposure.start_time = Time("2022-01-01T00:00:00")
    m.meta.exposure.end_time = Time("2022-01-01T01:00:00")
    model_path = _model_to_tmpfile(m, tmp_path)
    try:
        stp.add_wcs(model_path, tolerance=0, allow_default=True)
    except ValueError:
        pass  # This is what we want for the test.
    except Exception as e:
        pytest.skip(f"Live ENGDB service is not accessible.\nException={e}")

    with rdm.open(model_path) as result:
        assert result.meta.wcsinfo.ra_ref == result.meta.pointing.target_ra
        assert result.meta.wcsinfo.dec_ref == result.meta.pointing.target_dec


def test_change_engdb_url():
    """Test changing the engineering database by call for success.

    The given time and database should not find any values.
    """
    with pytest.raises(ValueError):
        stp.get_pointing(
            Time("2015-06-15"), Time("2015-06-17"), engdb_url=engdb_mast.MAST_BASE_URL
        )


def test_change_engdb_url_fail():
    """Test changing the engineering database by call"""
    with pytest.raises(ValueError):
        stp.get_pointing(
            Time(STARTTIME, format="isot"),
            Time(ENDTIME, format="isot"),
            engdb_url="https://nonexistent.fake.example",
        )


def test_get_pointing(engdb):
    """Ensure that the averaging works."""

    obstime, q = stp.get_pointing(STARTTIME, ENDTIME)

    assert np.isclose(obstime.value, OBSTIME_EXPECTED.value)
    assert np.allclose(q, Q_EXPECTED)


def test_get_pointing_fail(engdb):
    with pytest.raises(ValueError):
        obstime, q = stp.get_pointing(BADSTARTTIME, BADENDTIME)


def test_get_pointing_list(engdb):
    results = stp.get_pointing(
        STARTTIME.mjd, ENDTIME.mjd, reduce_func=stp.all_pointings
    )
    assert isinstance(results, list)
    assert len(results) > 0
    assert np.isclose(results[0].q, Q_EXPECTED).all()
    assert STARTTIME <= results[0].obstime <= ENDTIME


def test_logging(caplog, engdb):
    stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert "Determining pointing between observations times" in caplog.text
    assert "Telemetry search tolerance" in caplog.text
    assert "Reduction function" in caplog.text
    assert "Querying engineering DB" in caplog.text


@pytest.mark.parametrize("wcs_type", ["wcsinfo", "vinfo"])
def test_wcs(calc_wcs, wcs_type):
    """Ensure WCS information is correct

    Parameters
    ----------
    calc_wcs : pytest.fixture
        Equates to tuple (wcsinfo, vinfo, transforms, t_pars)

    wcs_type : str
        Which particular WCS information to test
    """
    wcsinfo, vinfo, _, t_pars = calc_wcs
    wcs = {"wcsinfo": wcsinfo, "vinfo": vinfo}
    with asdf.open(DATA_PATH / "wcs.asdf") as af:
        expected = af.tree[wcs_type]

    assert expected == wcs[wcs_type]._asdict()


def test_strict_pointing(science_raw_model, tmp_path):
    """Test failure on strict pointing"""
    science_raw_model.meta.exposure.end_time = STARTTIME
    with pytest.raises(ValueError):
        stp.add_wcs(_model_to_tmpfile(science_raw_model, tmp_path), tolerance=0)


@pytest.mark.parametrize(
    "matrix", [matrix for matrix in dataclasses.fields(stp.Transforms())]
)
def test_transforms(calc_wcs, matrix):
    """Ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare
    """
    _, _, transforms, t_pars = calc_wcs
    _test_transforms(transforms, t_pars, matrix.name)


def test_transform_serialize(calc_wcs, tmp_path):
    """Test serialization of Transforms"""
    _, _, transforms, _ = calc_wcs

    path = tmp_path / "transforms.asdf"
    transforms.write_to_asdf(path)
    from_asdf = stp.Transforms.from_asdf(path)

    assert isinstance(from_asdf, stp.Transforms)
    assert str(transforms) == str(from_asdf)


# ######################
# Utilities and fixtures
# ######################
def _make_t_pars(detector="WFI02"):
    """Setup initial Transforms Parameters

    This set is Visit 1 provided by T.Sohn in a demonstration notebook.

    Parameters
    ==========
    detector : str
        Detector in use.
    """
    t_pars = stp.TransformParameters()

    t_pars.detector = detector

    t_pars.guide_star_wcs = stp.WCSRef(ra=270.0, dec=66.5607, pa=None)
    t_pars.pointing = stp.Pointing(
        q=np.array([0.709433, -0.291814, 0.641319, 0.016107]),
    )

    return t_pars


def _model_to_tmpfile(m, tmp_path, fname="file.asdf"):
    """Save a DataModel to a general temp file"""
    file_path = tmp_path / fname
    m.save(file_path)
    m.close()
    return file_path


def _test_transforms(transforms, t_pars, matrix):
    """Private function to ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare
    """
    expected_tforms = stp.Transforms.from_asdf(DATA_PATH / "transforms.asdf")
    expected_value = getattr(expected_tforms, matrix)

    value = getattr(transforms, matrix)
    if expected_value is None:
        assert value is None
    else:
        assert np.allclose(value, expected_value)


@pytest.fixture(scope='module')
def engdb():
    """Setup the service to operate through the mock service"""
    try:
        engdb = engdb_tools.ENGDB_Service()
    except RuntimeError as exception:
        pytest.skip(f"Engineering database unvailable: {exception}.")
    yield engdb


@pytest.fixture
def calc_wcs(tmp_path_factory, engdb):
    """Calculate full transforms and WCS info"""
    t_pars = _make_t_pars()

    # Calculate the transforms and WCS information
    wcsinfo, vinfo, transforms = stp.calc_wcs(t_pars)

    # Save all for later examination.
    transforms_path = tmp_path_factory.mktemp("transforms")
    transforms.write_to_asdf(transforms_path / "transforms.asdf")
    wcs_asdf_file = asdf.AsdfFile(
        {"wcsinfo": wcsinfo._asdict(), "vinfo": vinfo._asdict()}
    )
    wcs_asdf_file.write_to(transforms_path / "wcs.asdf")

    return wcsinfo, vinfo, transforms, t_pars


@pytest.fixture
def science_raw_model():
    """Create the base model for testing"""
    m = rdm.datamodels.ScienceRawModel.create_fake_data(
        {
            "meta": {
                "exposure": {"start_time": STARTTIME, "end_time": ENDTIME},
                "pointing": {"target_ra": TARG_RA, "target_dec": TARG_DEC},
            }
        }
    )

    return m
