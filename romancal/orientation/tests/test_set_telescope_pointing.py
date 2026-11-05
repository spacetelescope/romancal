"""
Test suite for set_telescope_pointing
"""

import logging
from pathlib import Path

import asdf
import numpy as np
import pytest
import roman_datamodels as rdm
from astropy.time import Time

from romancal.lib.engdb import engdb_mast, engdb_tools
from romancal.orientation import _lib as olib
from romancal.orientation import _pointing as plib
from romancal.orientation import _transforms as tlib
from romancal.orientation import _wcs as wlib
from romancal.orientation import set_telescope_pointing as stp

# Ensure that `set_telescope_pointing` logs.
stp.logger.setLevel(logging.DEBUG)
stp.logger.addHandler(logging.StreamHandler())

# Get the mock databases
DATA_PATH = Path(__file__).parent / "data"

# Setup for some times
STARTTIME = Time("2026-09-22T06:18:21.121", format="isot")
ENDTIME = Time("2026-09-22T06:19:30.696", format="isot")
BADSTARTTIME = Time("2020-02-02T02:02:02", format="isot")
BADENDTIME = Time("2020-02-02T02:12:02", format="isot")
DEFAULT_RADECREF = (80.73771580167708, -69.37882680328173)

# Header defaults
TARG_RA = 270.0
TARG_DEC = 66.0

# Default set of transform parameters
# These are the set of initial conditions used in the notebook `demo_vacorrection` as
# sent by T.Sohn to J.Eisenhamer on 2015-11-19. The transforms produced are represented
# in the truth file `data/transforms.asdf` and the resulting WCS results in `data/wcs.asdf`
TRANSFORM_KWARGS = {
    "aperture": "WFI01_FULL",
    "bam_ref": DATA_PATH / "roman_wfi_bam_0007.asdf",
    "distortion_ref": DATA_PATH / "roman_wfi_distortion_0078.asdf",
    "pointing": plib.Pointing(
        obstime=Time(1790057901.121, format="unix"),
        q=np.array(
            [
                -0.7118218303725413,
                -0.1163575361747678,
                0.6786490188532825,
                0.13858612769623332,
            ]
        ),
    ),
    "target_aperture": "WFI_CEN",
    "velocity": (0.35121860071119626, 27.235876615565928, 11.79987217049065),
}

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
        assert np.allclose(
            result.meta.pointing.quaternion, TRANSFORM_KWARGS["pointing"].q
        )


def test_add_wcs_default_from_model(science_raw_model, tmp_path):
    """Handle when no pointing exists and the default is used."""
    m = science_raw_model
    m.meta.exposure.start_time = Time("2022-01-01T00:00:00")
    m.meta.exposure.end_time = Time("2022-01-01T01:00:00")
    m.meta.pointing.quaternion = tuple(TRANSFORM_KWARGS["pointing"].q)
    model_path = _model_to_tmpfile(m, tmp_path)
    stp.add_wcs(
        model_path,
        tolerance=0,
        allow_default=True,
    )

    with rdm.open(model_path) as result:
        assert np.isclose(result.meta.wcsinfo.ra_ref, DEFAULT_RADECREF[0])
        assert np.isclose(result.meta.wcsinfo.dec_ref, DEFAULT_RADECREF[1])
        assert np.allclose(
            result.meta.pointing.quaternion, TRANSFORM_KWARGS["pointing"].q
        )


def test_change_base_url():
    """Test changing the engineering database by call for success.

    The given time and database should not find any values.
    """
    service_kwargs = {"service": "mast", "eng_base_url": engdb_mast.MAST_BASE_URL}
    with pytest.raises(ValueError):
        plib.get_pointing(
            Time("2015-06-15"), Time("2015-06-17"), service_kwargs=service_kwargs
        )


def test_change_base_url_fail():
    """Test changing the engineering database by call"""
    service_kwargs = {
        "service": "mast",
        "eng_base_url": "https://nonexistent.fake.example",
    }
    with pytest.raises(ValueError):
        plib.get_pointing(
            Time(STARTTIME, format="isot"),
            Time(ENDTIME, format="isot"),
            service_kwargs=service_kwargs,
        )


@pytest.mark.parametrize('override, detector, target_aperture, partial_result',
                         [
                             # All references pulled from CRDS with nothing for WFI_CEN target.
                             ({'bam_ref': None, 'distortion_ref': None, 'distortion_target_ref': None},
                              'WFI01', 'WFI_CEN',
                              {'bam_ref': 'roman', 'distortion_ref': 'roman', 'distortion_target_ref': 'A'}),
                             # Overrides.
                             ({'bam_ref': 'bam_local', 'distortion_ref': 'dist_local', 'distortion_target_ref': 'targ_local'},
                              'WFI01', 'WFI_CEN',
                              {'bam_ref': 'bam', 'distortion_ref': 'dist', 'distortion_target_ref': 'A'}),
                             # Different target aperture.
                             ({'bam_ref': 'bam_local', 'distortion_ref': 'dist_local', 'distortion_target_ref': None},
                              'WFI01', 'WFI18_FULL',
                              {'bam_ref': 'bam', 'distortion_ref': 'dist', 'distortion_target_ref': 'roman'}),
                             # Different target aperture with override.
                             ({'bam_ref': 'bam_local', 'distortion_ref': 'dist_local', 'distortion_target_ref': 'targ_local'},
                              'WFI01', 'WFI18_FULL',
                              {'bam_ref': 'bam', 'distortion_ref': 'dist', 'distortion_target_ref': 'targ'}),
                         ])
def test_fetch_references(override, detector, target_aperture, partial_result):
    """Test fetching references."""
    t_pars = tlib.TransformParameters(**TRANSFORM_KWARGS)
    t_pars.target_aperture = target_aperture
    for ref in override:
        setattr(t_pars, ref, override[ref])
    t_pars.crds_pars = olib.mk_crds_pars(t_pars.pointing.obstime, detector=detector)

    olib.fetch_references(t_pars)

    for ref in partial_result:
        assert Path(getattr(t_pars, ref)).name.startswith(partial_result[ref])


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_get_mnemonics():
    """Test getting mnemonics"""
    try:
        mnemonics = plib.get_mnemonics(STARTTIME, ENDTIME, 60)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))

    assert len(mnemonics) == len(plib.COARSE_MNEMONICS)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_get_pointing():
    """Ensure that the averaging works."""
    try:
        pointing = plib.get_pointing(STARTTIME, ENDTIME)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))

    assert np.isclose(pointing.obstime.value, 1790057935.6530108)
    assert np.allclose(
        pointing.q, np.array([-0.71182183, -0.11635754, 0.67864902, 0.13858613])
    )


def test_get_pointing_fail():
    with pytest.raises(ValueError):
        obstime, q = plib.get_pointing(BADSTARTTIME, BADENDTIME)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_get_pointing_list():
    """Test pointing collection

    Note: The expected quaternion is from mastdev during Q26B21 development.
    This will most likely change again.
    """
    try:
        results = plib.get_pointing(STARTTIME, ENDTIME, reduce_func=plib.all_pointings)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))
    assert isinstance(results, list)
    assert len(results) > 0
    assert np.isclose(
        results[0].q,
        np.array([-0.71182157, -0.11635775, 0.67864924, 0.13858621]),
        rtol=1.0e-2,
    ).all()
    assert STARTTIME <= results[0].obstime <= ENDTIME


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_mnemonics_chronologically():
    """Test ordering mnemonics chronologically"""
    try:
        mnemonics = plib.get_mnemonics(STARTTIME, ENDTIME, 60)
    except ValueError as exception:
        pytest.xfail(reason=str(exception))
    ordered = plib.mnemonics_chronologically(mnemonics)

    assert len(ordered) > 1

    first = ordered[0]
    assert isinstance(first[0], Time)
    assert isinstance(first[1], dict)
    assert len(first[1]) >= len(plib.COARSE_MNEMONICS)


@pytest.mark.skipif(NO_ENGDB, reason="No engineering database available")
def test_logging(caplog):
    try:
        plib.get_pointing(STARTTIME, ENDTIME)
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
        )
    )

    assert expected == set(plib.COARSE_MNEMONICS)


def test_wcs(calc_wcs):
    """Ensure WCS information is correct

    Parameters
    ----------
    calc_wcs : pytest.fixture
        Equates to tuple (wcsinfo, vinfo, transforms, t_pars)
    """
    wcs_meta, t_pars = calc_wcs
    with asdf.open(DATA_PATH / "wcs.asdf") as af:
        expected = af.tree["wcs_meta"]

    wcs_meta_flat = _flatten_dict(wcs_meta)
    expected_flat = _flatten_dict(expected)
    for key in wcs_meta_flat:
        try:
            float(wcs_meta_flat[key])
        except (TypeError, ValueError):
            assert expected_flat[key] == wcs_meta_flat[key]
        else:
            assert np.isclose(expected_flat[key], wcs_meta_flat[key]), (
                f"Key {key} differs expected {expected_flat[key]} calculated {wcs_meta_flat[key]}"
            )


def test_strict_pointing(science_raw_model, tmp_path):
    """Test failure on strict pointing"""
    science_raw_model.meta.exposure.end_time = STARTTIME
    with pytest.raises(ValueError):
        stp.add_wcs(_model_to_tmpfile(science_raw_model, tmp_path), tolerance=0)


@pytest.mark.parametrize(
    "attr, expected",
    [
        ("meta.pointing.pa_aperture", 37.24854770366899),
        ("meta.pointing.pointing_engineering_source", "CALCULATED"),
        ("meta.pointing.ra_v1", 81.6686302053683),
        ("meta.pointing.dec_v1", -69.0573058072646),
        ("meta.pointing.pa_v3", 96.49949482075951),
        ("meta.pointing.target_aperture", "WFI_CEN"),
        ("meta.pointing.target_ra", 80.71805432662013),
        ("meta.pointing.target_dec", -69.48828716074138),
        (
            "meta.pointing.quaternion",
            [
                -0.7118218303725413,
                -0.1163575361747678,
                0.6786490188532825,
                0.13858612769623332,
            ],
        ),
        (
            "meta.pointing.bam_quaternion",
            [
                -0.18541749359303772,
                0.6841636922799866,
                -0.17893482160935997,
                0.6822922576514541,
            ],
        ),
        ("meta.wcsinfo.aperture_name", "WFI01_FULL"),
        ("meta.wcsinfo.v2_ref", 1454.084),
        ("meta.wcsinfo.v3_ref", -1163.304),
        ("meta.wcsinfo.vparity", -1),
        ("meta.wcsinfo.v3yangle", -59.99960999999999),
        ("meta.wcsinfo.ra_ref", 80.62520113957764),
        ("meta.wcsinfo.dec_ref", -69.41891635642607),
        ("meta.wcsinfo.roll_ref", 97.47628097051559),
        (
            "meta.wcsinfo.s_region",
            "POLYGON ICRS 80.660744006 -69.505616219 80.379384214 -69.429205499 80.590434179 -69.331579227 80.872932680 -69.407979480",
        ),
    ],
)
def test_update_meta(attr, expected, updated_model):
    """Ensure that all the expected meta are updated"""

    value = _getattrpath(updated_model, attr)
    if isinstance(value, (list, tuple, np.ndarray)):
        assert np.allclose(value, expected)
    elif isinstance(value, float):
        assert np.isclose(value, expected)
    else:
        assert value == expected


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
    wcs_meta = wlib.calc_wcs(t_pars)

    # Save all for later examination.
    wcs_meta_path = tmp_path_factory.mktemp("transforms")
    wcs_asdf_file = asdf.AsdfFile({"wcs_meta": wcs_meta})
    wcs_asdf_file.write_to(wcs_meta_path / "wcs.asdf")

    return wcs_meta, t_pars


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


@pytest.fixture(scope="module")
def updated_model(calc_wcs):
    """Update a model with a orientation calculation"""
    wcs_meta, t_pars = calc_wcs

    m = rdm.datamodels.ScienceRawModel.create_fake_data(
        {
            "meta": {
                "ephemeris": {
                    "velocity_x": TRANSFORM_KWARGS["velocity"][0],
                    "velocity_y": TRANSFORM_KWARGS["velocity"][1],
                    "velocity_z": TRANSFORM_KWARGS["velocity"][2],
                },
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
    olib.update_meta(m, t_pars, wcs_meta, "CALCULATED")
    return m


def _flatten_dict(d, parent_key="", sep="."):
    """Flattens a nested dictionary, joining keys with a separator."""
    items = {}
    for key, value in d.items():
        # Build the new key path
        new_key = f"{parent_key}{sep}{key}" if parent_key else key
        if isinstance(value, dict):
            items.update(_flatten_dict(value, new_key, sep=sep))
        else:
            items[new_key] = value
    return items


def _getattrpath(obj, path, default=None):
    """Retrieve an attribute from its fully qualified path"""
    attribute, _, rest = path.partition(".")
    value = getattr(obj, attribute, default)
    if value is not None and rest:
        value = _getattrpath(value, rest, default)
    return value


def _make_t_pars(**transform_kwargs):
    """Setup initial Transforms Parameters

    Parameters
    ==========
    transform_kwargs : dict
        dict to use to initialize the `TransformParameters` object.
        See `TransformParameters` for more information.`
    """
    t_pars = tlib.TransformParameters(**transform_kwargs)

    # Force MAST service
    t_pars.service_kwargs = {"service": "mast"}

    return t_pars


def _model_to_tmpfile(m, tmp_path, fname="file.asdf"):
    """Save a DataModel to a general temp file"""
    file_path = tmp_path / fname
    m.save(file_path)
    m.close()
    return file_path


def _test_transforms(transforms, t_pars, matrix, expected_fname="transforms.asdf"):
    """Private function to ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare
    """
    expected_tforms = tlib.Transforms.from_asdf(DATA_PATH / expected_fname)
    expected_value = getattr(expected_tforms, matrix)

    value = getattr(transforms, matrix)
    if expected_value is None:
        assert value is None
    else:
        assert np.allclose(value, expected_value)
