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
    t_pars = make_t_pars()
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
def make_t_pars(detector='WFI02'):
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
    t_pars = make_t_pars()

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
def data_file(tmp_path):
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '2017-01-01'
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x = -25.021
    model.meta.ephemeris.velocity_y = -16.507
    model.meta.ephemeris.velocity_z = -7.187

    file_path = tmp_path / 'file.fits'
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture(params=[('good_model', True), ('bad_model', False), ('fits_nomodel', False)])
def file_case(request, tmp_path):
    """Generate files with different model states"""
    case, allow = request.param

    if case == 'good_model':
        # Make a model that will always succeed
        model = datamodels.Level1bModel((10, 10, 10, 10))
        path = tmp_path / 'level1bmodel.fits'
        model.save(path)
    elif case == 'bad_model':
        # Make a model that will fail if not allowed
        model = datamodels.IFUCubeModel((10, 10, 10))
        path = tmp_path / 'image.fits'
        model.save(path)
    elif case == 'fits_nomodel':
        # Create just a plain anything FITS
        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        path = tmp_path / 'empty.fits'
        hdul.writeto(path)
    else:
        assert False, f'Cannot produce a file for {case}'

    return path, allow


@pytest.fixture
def data_file_nosiaf(tmp_path):
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.aperture.name = "UNKNOWN"
    model.meta.observation.date = '2017-01-01'

    file_path = tmp_path / 'fits_nosiaf.fits'
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_fromsim(tmp_path):
    """Create data using times that were executed during a simulation using the OTB Simulator"""
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = Time('2022-02-02T22:24:58.942').mjd
    model.meta.exposure.end_time = Time('2022-02-02T22:26:24.836').mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.guidestar.gs_pcs_mode = 'COARSE'
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '2017-01-01'
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x_bary = -25.021
    model.meta.ephemeris.velocity_y_bary = -16.507
    model.meta.ephemeris.velocity_z_bary = -7.187

    file_path = tmp_path / 'file_fromsim.fits'
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_acq1(tmp_path):
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = STARTTIME.mjd
    model.meta.exposure.end_time = ENDTIME.mjd
    model.meta.target.ra = TARG_RA
    model.meta.target.dec = TARG_DEC
    model.meta.guidestar.gs_ra = TARG_RA + 0.0001
    model.meta.guidestar.gs_dec = TARG_DEC + 0.0001
    model.meta.aperture.name = "FGS2_FULL"
    model.meta.observation.date = '2017-01-01'
    model.meta.exposure.type = "FGS_ACQ1"
    model.meta.ephemeris.velocity_x = -25.021
    model.meta.ephemeris.velocity_y = -16.507
    model.meta.ephemeris.velocity_z = -7.187

    file_path = tmp_path / 'file.fits'
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def data_file_moving_target(tmp_path):
    """Example data from simulation."""
    # Values are from simulated data file jw00634_nrcblong_mttest_uncal.fits
    model = datamodels.Level1bModel()
    model.meta.exposure.start_time = 58738.82598848102
    model.meta.exposure.end_time = 58738.82747969907
    model.meta.exposure.mid_time = 58738.82673409005
    model.meta.target.ra = 0.0
    model.meta.target.dec = 0.0
    model.meta.guidestar.gs_ra = 0.0001
    model.meta.guidestar.gs_dec = 0.0001
    model.meta.aperture.name = "MIRIM_FULL"
    model.meta.observation.date = '2019-09-12'
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.ephemeris.velocity_x = 0.00651191175424979
    model.meta.ephemeris.velocity_y = 0.160769793796114
    model.meta.ephemeris.velocity_z = 0.147663026601154

    model.meta.target.type = 'MOVING'
    model.meta.moving_target = None

    times = ['2019-09-12T19:49:25.405', '2019-09-12T19:50:29.825', '2019-09-12T19:51:34.246']
    apparent_ra = [0.0, 6.2e-5, 1.24e-4]
    apparent_dec = [-6.2e-5,  0.0,  3.0e-5]
    default = [0.0, 0.0, 0.0]
    col_names = [item['name'] for item in model.schema['properties']['moving_target']['datatype']]
    mt_table = Table([times, apparent_ra, apparent_dec],
                     names=('time', 'mt_apparent_RA', 'mt_apparent_Dec'))
    for column in col_names:
        if column not in {'time', 'mt_apparent_RA', 'mt_apparent_Dec'}:
            mt_table.add_column(default, name=column)
    model.moving_target = mt_table.as_array()

    file_path = tmp_path / 'file.fits'
    model.save(file_path)
    model.close()
    yield file_path


@pytest.fixture
def science_raw_model():
    """Create the base model for testing"""
    m = rdm.datamodels.ScienceRawModel.create_fake_data({'meta': {
        'exposure': {'start_time': STARTTIME, 'end_time': ENDTIME},
        'pointing': {'target_ra': TARG_RA, 'target_dec': TARG_DEC},
    }})

    return m
