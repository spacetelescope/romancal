"""
Test suite for set_telescope_pointing
"""

import dataclasses
import logging
from pathlib import Path

import asdf
import numpy as np
import pytest
import roman_datamodels as rdm
from astropy.time import Time

from romancal.lib.engdb import engdb_mast, engdb_tools
from romancal.orientation import set_telescope_pointing as stp

# pysiaf is not a required dependency. If not present, ignore all this.
pysiaf = pytest.importorskip("pysiaf")

# Ensure that `set_telescope_pointing` logs.
stp.logger.setLevel(logging.DEBUG)
stp.logger.addHandler(logging.StreamHandler())

# Setup for some times
STARTTIME = Time("2027-03-11T13:26:56.571", format="isot")
ENDTIME = Time("2027-03-11T13:27:56.658", format="isot")
BADSTARTTIME = Time("2020-02-02T02:02:02", format="isot")
BADENDTIME = Time("2020-02-02T02:12:02", format="isot")
DEFAULT_RADECREF = (149.9709175757851, 86.7401422250309)

# Header defaults
TARG_RA = 270.0
TARG_DEC = 66.01

# Default set of transform parameters
TRANSFORM_KWARGS = {
    "aperture": "WFI01_FULL",
    "gscommanded": (916.4728835141, -186.8939737044),
    "pointing": stp.Pointing(
        fgs_q=np.array(
            [
                -0.18596734175399293,
                0.6837984564491885,
                -0.1800546332580956,
                0.6822141509826322,
            ]
        ),
        obstime=Time(1804771646.4242268, format="unix"),
        q=np.array([-0.33879082, 0.62326573, -0.36611627, 0.60226181]),
    ),
    "velocity": (-5.473753741352352, -27.480586797035414, -11.875972151015253),
}

# Get the mock databases
DATA_PATH = Path(__file__).parent / "data"

# Meta attributes for test comparisons
METAS_EQUALITY = [
    "meta.pointing.pointing_engineering_source",
]
METAS_ISCLOSE = [
    "meta.wcsinfo.dec_ref",
    "meta.wcsinfo.ra_ref",
    "meta.wcsinfo.roll_ref",
    "meta.wcsinfo.v2_ref",
    "meta.wcsinfo.v3_ref",
    "meta.pointing.ra_v1",
    "meta.pointing.dec_v1",
    "meta.pointing.pa_v3",
]

# Some tests depend on the default engineering database to be available. Check now.
NO_ENGDB = False
try:
    engdb_tools.engdb_service()
except engdb_tools.EXPECTED_ERRORS:
    NO_ENGDB = True


def test_add_wcs_default(science_raw_model, tmp_path):
    """Handle when no pointing exists and the default is used."""
    m = science_raw_model
    m.meta.exposure.start_time = Time("2022-01-01T00:00:00")
    m.meta.exposure.end_time = Time("2022-01-01T01:00:00")
    model_path = _model_to_tmpfile(m, tmp_path)
    stp.add_wcs(
        model_path,
        tolerance=0,
        allow_default=True,
        default_quaternion=TRANSFORM_KWARGS["pointing"].q,
    )

    with rdm.open(model_path) as result:
        assert np.isclose(result.meta.wcsinfo.ra_ref, DEFAULT_RADECREF[0])
        assert np.isclose(result.meta.wcsinfo.dec_ref, DEFAULT_RADECREF[1])


def test_change_base_url():
    """Test changing the engineering database by call for success.

    The given time and database should not find any values.
    """
    service_kwargs = {"service": "mast", "eng_base_url": engdb_mast.MAST_BASE_URL}
    with pytest.raises(ValueError):
        stp.get_pointing(
            Time("2015-06-15"), Time("2015-06-17"), service_kwargs=service_kwargs
        )


def test_change_base_url_fail():
    """Test changing the engineering database by call"""
    service_kwargs = {
        "service": "mast",
        "eng_base_url": "https://nonexistent.fake.example",
    }
    with pytest.raises(ValueError):
        stp.get_pointing(
            Time(STARTTIME, format="isot"),
            Time(ENDTIME, format="isot"),
            service_kwargs=service_kwargs,
        )


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_get_mnemonics():
    """Test getting mnemonics"""
    try:
        mnemonics = stp.get_mnemonics(STARTTIME, ENDTIME, 60)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))

    assert len(mnemonics) == len(stp.COARSE_MNEMONICS)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_get_pointing():
    """Ensure that the averaging works.

    Note: The expected quaternion is from mastdev during Q26B20 development.
    This will most likely change again.
    """
    try:
        pointing = stp.get_pointing(STARTTIME, ENDTIME)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))

    assert np.isclose(
        pointing.obstime.value, TRANSFORM_KWARGS["pointing"].obstime.value
    )
    assert np.allclose(pointing.q, TRANSFORM_KWARGS["pointing"].q)


def test_get_pointing_fail():
    with pytest.raises(ValueError):
        obstime, q = stp.get_pointing(BADSTARTTIME, BADENDTIME)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_get_pointing_list():
    """Test pointing collection

    Note: The expected quaternion is from mastdev during Q26B20 development.
    This will most likely change again.
    """
    try:
        results = stp.get_pointing(STARTTIME, ENDTIME, reduce_func=stp.all_pointings)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))
    assert isinstance(results, list)
    assert len(results) > 0
    assert np.isclose(results[0].q, TRANSFORM_KWARGS["pointing"].q, rtol=1.0e-2).all()
    assert STARTTIME <= results[0].obstime <= ENDTIME


def test_hv_to_fgs():
    """Test conversion from HV frame to FGS frame"""
    hv = TRANSFORM_KWARGS["gscommanded"]
    fgs_expected = (346.06680318732685, -148.7528870949794)

    fgs = stp.hv_to_fgs("WFI01_FULL", *hv, pysiaf)

    assert np.allclose(fgs, fgs_expected)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_mnemonics_chronologically():
    """Test ordering mnemonics chronologically"""
    try:
        mnemonics = stp.get_mnemonics(STARTTIME, ENDTIME, 60)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))
    ordered = stp.mnemonics_chronologically(mnemonics)

    assert len(ordered) > 1

    first = ordered[0]
    assert isinstance(first[0], Time)
    assert isinstance(first[1], dict)
    assert len(first[1]) >= len(stp.COARSE_MNEMONICS_QUATERNION_ECI)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_logging(caplog):
    try:
        stp.get_pointing(STARTTIME, ENDTIME)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))
    assert "Determining pointing between observations times" in caplog.text
    assert "Telemetry search tolerance" in caplog.text
    assert "Reduction function" in caplog.text
    assert "Querying engineering DB" in caplog.text


def test_mnemonic_list():
    """Ensure the mnemonic list is as expected"""
    expected = set(
        (
            "SCF_AC_SDR_QBJ_1",
            "SCF_AC_SDR_QBJ_2",
            "SCF_AC_SDR_QBJ_3",
            "SCF_AC_SDR_QBJ_4",
            "SCF_AC_EST_FGS_qbr1",
            "SCF_AC_EST_FGS_qbr2",
            "SCF_AC_EST_FGS_qbr3",
            "SCF_AC_EST_FGS_qbr4",
            "SCF_AC_FGS_TBL_Qb1",
            "SCF_AC_FGS_TBL_Qb2",
            "SCF_AC_FGS_TBL_Qb3",
            "SCF_AC_FGS_TBL_Qb4",
        )
    )

    assert expected == set(stp.COARSE_MNEMONICS)


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

    wcs_dict = wcs[wcs_type]._asdict()
    for key in wcs_dict:
        if key != "s_region":
            assert np.isclose(expected[key], wcs_dict[key]), (
                f"Key {key} differs expected {expected[key]} calculated {wcs_dict[key]}"
            )
        else:
            assert expected[key] == wcs_dict[key]


def test_strict_pointing(science_raw_model, tmp_path):
    """Test failure on strict pointing"""
    science_raw_model.meta.exposure.end_time = STARTTIME
    with pytest.raises(ValueError):
        stp.add_wcs(_model_to_tmpfile(science_raw_model, tmp_path), tolerance=0)


@pytest.mark.parametrize(
    "matrix",
    [matrix for matrix in dataclasses.fields(stp.Transforms())],
    ids=[matrix.name for matrix in dataclasses.fields(stp.Transforms())],
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
@pytest.fixture(scope="module")
def calc_wcs(tmp_path_factory):
    """Calculate full transforms and WCS info

    ***DEFINE WHERE DEFAULTS HAVE COME FROM***

    """
    t_pars = _make_t_pars(**TRANSFORM_KWARGS)

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
                "pointing": {
                    "target_aperture": "WFI_CEN",
                    "target_ra": TARG_RA,
                    "target_dec": TARG_DEC,
                },
                "wcsinfo": {"aperture_name": "WFI02_FULL"},
            }
        }
    )

    return m


def _make_t_pars(**transform_kwargs):
    """Setup initial Transforms Parameters

    Parameters
    ==========
    transform_kwargs : dict
        dict to use to initialize the `TransformParameters` object.
        See `TransformParameters` for more information.`
    """
    t_pars = stp.TransformParameters(**transform_kwargs)

    # Force MAST service
    t_pars.service_kwargs = {"service": "mast"}

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
