"""
Test suite for set_telescope_pointing
"""
import logging
import numpy as np
from pathlib import Path
import pytest
import warnings

pytest.importorskip('pysiaf')

import asdf

from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.time import Time  # noqa: E402

from romancal.lib import engdb_mast  # noqa: E402
from romancal.orientation import set_telescope_pointing as stp  # noqa: E402
from romancal.lib import siafdb  # noqa: E402
import roman_datamodels as rdm

# Ensure that `set_telescope_pointing` logs.
stp.logger.setLevel(logging.DEBUG)
stp.logger.addHandler(logging.StreamHandler())

# Setup mock engineering service
STARTTIME = Time('2027-03-23T19:20:40', format='isot')
ENDTIME = Time('2027-03-23T19:21:36', format='isot')

# Header defaults
TARG_RA = 270.0
TARG_DEC = 66.01

# Get the mock databases
DATA_PATH = Path(__file__).parent / 'data'

Q_EXPECTED = np.array([-0.69018802,  0.12195182, -0.695103  ,  0.15999998])
OBSTIME_EXPECTED = Time(1805829668.0276294, format='unix')

# Meta attributes for test comparisons
METAS_EQUALITY = ['meta.exposure.engineering_quality',
                  ]
METAS_ISCLOSE = [
    'meta.wcsinfo.roll_ref',
    'meta.wcsinfo.v2_ref',
    'meta.wcsinfo.v3_ref',
    'meta.wcsinfo.ra_ref',
    'meta.wcsinfo.dec_ref',
    'meta.pointing.ra_v1',
    'meta.pointing.dec_v1',
    'meta.pointing.pa_v3',
]


def test_add_wcs_default(science_raw_model, tmp_path):
    """Handle when no pointing exists and the default is used."""
    m = science_raw_model
    m.meta.exposure.start_time = Time('2022-01-01T00:00:00')
    m.meta.exposure.end_time = Time('2022-01-01T01:00:00')
    model_path = _model_to_tmpfile(m, tmp_path)
    try:
        stp.add_wcs(
            model_path, tolerance=0, allow_default=True
        )
    except ValueError:
        pass  # This is what we want for the test.
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            f'\nException={e}'
        )

    with rdm.open(model_path) as result:
        assert result.meta.wcsinfo.ra_ref == result.meta.pointing.target_ra
        assert result.meta.wcsinfo.dec_ref == result.meta.pointing.target_dec


def test_change_engdb_url():
    """Test changing the engineering database by call for success.

    The given time and database should not find any values.
    """
    with pytest.raises(ValueError):
        stp.get_pointing(
            Time('2015-06-15'),
            Time('2015-06-17'),
            engdb_url=engdb_mast.MAST_BASE_URL
        )


def test_change_engdb_url_fail():
    """Test changing the engineering database by call"""
    with pytest.raises(Exception):
        stp.get_pointing(
            Time(STARTTIME, format='isot'),
            Time(ENDTIME, format='isot'),
            engdb_url='http://nonexistent.fake.example'
        )


def test_get_pointing():
    """Ensure that the averaging works."""

    obstime, q = stp.get_pointing(
         STARTTIME.mjd,
         ENDTIME.mjd
    )

    assert np.isclose(obstime.value, OBSTIME_EXPECTED.value)
    assert np.allclose(q, Q_EXPECTED)


def test_get_pointing_fail():
    with pytest.raises(Exception):
        obstime, q  = stp.get_pointing(47892.0, 48256.0)


def test_get_pointing_list():
    results = stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd, reduce_func=stp.all_pointings)
    assert isinstance(results, list)
    assert len(results) > 0
    assert np.isclose(results[0].q, Q_EXPECTED).all()
    assert STARTTIME <= results[0].obstime <= ENDTIME


def test_logging(caplog):
    stp.get_pointing(STARTTIME.mjd, ENDTIME.mjd)
    assert 'Determining pointing between observations times' in caplog.text
    assert 'Telemetry search tolerance' in caplog.text
    assert 'Reduction function' in caplog.text
    assert 'Querying engineering DB' in caplog.text


@pytest.mark.parametrize(
    'method',
    [method for method in stp.Methods]
)
def test_method_string(method):
    """Ensure that the value of the method is the string representation"""
    assert f'{method}' == method.value


@pytest.mark.parametrize('wcs_type', ['wcsinfo', 'vinfo'])
def test_method_wcs(calc_method, wcs_type, truth_ext=''):
    """Ensure WCS information is correct

    Parameters
    ----------
    calc_method : pytest.fixture
        Equates to tuple (wcsinfo, vinfo, transforms, t_pars)

    wcs_type : str
        Which particular WCS information to test

    truth_ext : str
        Arbitrary extension to add to the truth file name.
    """
    wcsinfo, vinfo, _, t_pars = calc_method
    wcs = {'wcsinfo': wcsinfo, 'vinfo': vinfo}
    with asdf.open(DATA_PATH / f'wcs_{t_pars.method}{truth_ext}.asdf') as af:
        expected = af.tree[wcs_type]

    assert expected == wcs[wcs_type]._asdict()


@pytest.mark.parametrize(
    'attribute, expected',
    [('m_eci2b', 'overridden'), ('m_eci2v', 'untouched')]
)
def test_override(attribute, expected):
    """Test overriding of Transforms attributes"""
    overrides = stp.Transforms(m_eci2b='overridden')
    to_override = stp.Transforms(m_eci2b='original', m_eci2v='untouched', override=overrides)

    assert getattr(to_override, attribute) == expected


def test_override_calc_wcs():
    """Test matrix override in the full calculation"""
    t_pars = _make_t_pars()
    wcsinfo, vinfo, transforms = stp.calc_wcs(t_pars)

    override = stp.Transforms(m_eci2b=np.array([[-0.17690118, -0.39338551,  0.91934622],
                                                [-0.43470441, -0.82917048, -0.35143805],
                                                [ 0.90054523, -0.3971454 , 0.00710906]]))
    t_pars.override_transforms = override
    wcsinfo_new, vinfo_new, transforms_new = stp.calc_wcs(t_pars)

    assert vinfo_new != vinfo
    assert all(np.isclose(vinfo_new,
                          stp.WCSRef(ra=245.78706748976023, dec=66.83068216214627, pa=89.45804357482956)))


def test_strict_pointing(science_raw_model, tmp_path):
    """Test failure on strict pointing"""
    science_raw_model.meta.exposure.end_time = STARTTIME
    with pytest.raises(ValueError):
        stp.add_wcs(_model_to_tmpfile(science_raw_model, tmp_path), tolerance=0)


@pytest.mark.parametrize('matrix', [matrix for matrix in stp.Transforms()._fields])
def test_transforms(calc_method, matrix):
    """Ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare
    """
    _, _, transforms, t_pars = calc_method
    _test_transforms(transforms, t_pars, matrix)


def test_transform_serialize(calc_method, tmp_path):
    """Test serialization of Transforms"""
    _, _, transforms, _ = calc_method

    path = tmp_path / 'transforms.asdf'
    transforms.write_to_asdf(path)
    from_asdf = stp.Transforms.from_asdf(path)

    assert isinstance(from_asdf, stp.Transforms)
    assert str(transforms) == str(from_asdf)


# ######################
# Utilities and fixtures
# ######################
def _make_t_pars(detector='WFI02'):
    """Setup initial Transforms Parameters

    This set is Visit 1 provided by T.Sohn in a demonstration notebook.

    Parameters
    ==========
    detector : str
        Detector in use.
    """
    t_pars = stp.TransformParameters()

    t_pars.detector = detector

    t_pars.guide_star_wcs = stp.WCSRef(ra=270., dec=66.5607, pa=None)
    t_pars.pointing = stp.Pointing(
        q=np.array([0.709433, -0.291814,  0.641319, 0.016107]),
    )

    return t_pars


def _calc_coarse_202111_fgsid_idfunc(value):
    """Created test IDS for calc_coarse_202111_fgsid"""
    detector, fgsid_user, fgs_expected = value
    return f'{detector}-{fgsid_user}'


def _model_to_tmpfile(m, tmp_path, fname='file.asdf'):
    """Save a DataModel to a general temp file"""
    file_path = tmp_path / fname
    m.save(file_path)
    m.close()
    return file_path


def _test_transforms(transforms, t_pars, matrix, truth_ext=''):
    """Private function to ensure expected calculate of the specified matrix

    Parameters
    ----------
    transforms, t_pars : Transforms, TransformParameters
        The transforms and the parameters used to generate the transforms

    matrix : str
        The matrix to compare

    truth_ext : str
        Arbitrary extension to add to the truth file name.
    """
    expected_tforms = stp.Transforms.from_asdf(DATA_PATH / f'transforms_{t_pars.method}{truth_ext}.asdf')
    expected_value = getattr(expected_tforms, matrix)

    value = getattr(transforms, matrix)
    if expected_value is None:
        assert value is None
    else:
        assert np.allclose(value, expected_value)


@pytest.fixture(scope='module',
                params=[method for method in stp.Methods])
def calc_method(request, tmp_path_factory):
    """Calculate full transforms and WCS info for a Method
    """
    t_pars = _make_t_pars()

    # Set the method
    t_pars.method = request.param

    # Calculate the transforms and WCS information
    wcsinfo, vinfo, transforms = stp.calc_wcs(t_pars)

    # Save all for later examination.
    transforms_path = tmp_path_factory.mktemp('transforms')
    transforms.write_to_asdf(transforms_path / f'transforms_{request.param}.asdf')
    wcs_asdf_file = asdf.AsdfFile({'wcsinfo': wcsinfo._asdict(), 'vinfo': vinfo._asdict()})
    wcs_asdf_file.write_to(transforms_path / f'wcs_{request.param}.asdf')

    return wcsinfo, vinfo, transforms, t_pars


@pytest.fixture
def science_raw_model():
    """Create the base model for testing"""
    m = rdm.datamodels.ScienceRawModel.create_fake_data({'meta': {
        'exposure': {'start_time': STARTTIME, 'end_time': ENDTIME},
        'pointing': {'target_ra': TARG_RA, 'target_dec': TARG_DEC},
    }})

    return m
