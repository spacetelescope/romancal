import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.maker_utils import mk_level2_image

from romancal.lib import dqflags
from romancal.source_catalog.source_catalog_step import SourceCatalogStep

DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
SATURATED = dqflags.pixel["SATURATED"]


def mk_image_model(
    rate_mean=0,
    rate_std=1e-4,
    image_shape=(100, 100),
    rng=np.random.default_rng(619),
):
    l2 = mk_level2_image(shape=image_shape)
    l2_im = ImageModel(l2)
    l2_im.data = u.Quantity(
        rng.normal(loc=rate_mean, scale=rate_std, size=l2_im.data.shape).astype(
            np.float32
        ),
        l2_im.data.unit,
    )

    # fake a background until `rad` implements the schema:
    l2_im.meta["background"] = dict(level=None, subtracted=False, method=None)

    l2_im.meta.cal_step["skymatch"] = "INCOMPLETE"
    return l2_im


@pytest.fixture
def wfi_model():
    units = u.electron / u.s
    rng = np.random.default_rng(seed=123)
    data = rng.normal(0, 0.5, size=(101, 101))
    data[20:80, 10:20] = 1.4
    data[20:30, 20:45] = 1.4
    data[20:80, 55:65] = 7.2
    data[70:80, 65:87] = 7.2
    data[45:55, 65:87] = 7.2
    data[20:30, 65:87] = 7.2
    data[55:75, 82:92] = 7.2
    data[25:45, 82:92] = 7.2

    err = np.abs(data) / 10.0

    model = mk_image_model(0, 0.5, image_shape=(101, 101), rng=rng)
    model.data = data << units
    model.err = err << units

    # model.meta.bunit_data = 'e/s'
    # model.meta.bunit_err = 'e/s'

    return model


# TODO (@bmorris3): replace or remove
# @pytest.fixture
# def wfi_model_without_apcorr():
#     rng = np.random.default_rng(seed=123)
#     data = rng.normal(0, 0.5, size=(101, 101))
#     data[20:80, 10:20] = 1.4
#     data[20:30, 20:45] = 1.4
#     data[20:80, 55:65] = 7.2
#     data[70:80, 65:87] = 7.2
#     data[45:55, 65:87] = 7.2
#     data[20:30, 65:87] = 7.2
#     data[55:75, 82:92] = 7.2
#     data[25:45, 82:92] = 7.2
#
#     wht = np.ones(data.shape)
#     wht[0:10, :] = 0.
#     err = np.abs(data) / 10.
#     model = ImageModel(data, wht=wht, err=err)
#     model.meta.bunit_data = 'MJy/sr'
#     model.meta.bunit_err = 'MJy/sr'
#     model.meta.photometry.pixelarea_steradians = 1.0
#     model.meta.wcs = make_gwcs(data.shape)
#     model.meta.wcsinfo = {
#         'ctype1': 'RA---TAN',
#         'ctype2': 'DEC--TAN',
#         'dec_ref': 11.99875540218638,
#         'ra_ref': 22.02351763251896,
#         'roll_ref': 0.005076934167039675,
#         'v2_ref': 86.039011,
#         'v3_ref': -493.385704,
#         'v3yangle': -0.07385127,
#         'vparity': -1,
#         'wcsaxes': 2,
#         'crpix1': 50,
#         'crpix2': 50}
#     model.meta.instrument = {
#         'channel': 'LONG',
#         'detector': 'NRCALONG',
#         'filter': 'F2550WR',
#         'lamp_mode': 'NONE',
#         'module': 'A',
#         'name': 'NIRCAM',
#         'pupil': 'CLEAR'}
#     model.meta.exposure.type = 'NRC_IMAGE'
#     model.meta.observation.date = '2021-01-01'
#     model.meta.observation.time = '00:00:00'
#
#     return model


@pytest.mark.parametrize("npixels, nsources", ((5, 2), (1000, 1), (5000, 0)))
def test_source_catalog(wfi_model, npixels, nsources):

    step = SourceCatalogStep(
        snr_threshold=0.5,
        npixels=npixels,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
    )
    cat = step.run(wfi_model)
    if cat is None:
        assert nsources == 0
    else:
        assert len(cat) == nsources
        min_snr = np.min(cat["isophotal_flux"] / cat["isophotal_flux_err"])
        assert min_snr >= 100


# TODO (@bmorris3): replace or remove
# def test_model_without_apcorr_data(wfi_model_without_apcorr):
#     step = SourceCatalogStep(save_results=False)
#     cat = step.run(wfi_model_without_apcorr)
#     assert cat is None


def test_input_model_reset(wfi_model):
    """Changes to input model data are made in SourceCatalogStep - make sure that
    these changes are correctly reversed so the input model data/err arrays
    remain unchanged after processing (to avoid copying datamodel),
    and that the units are in MJy/sr before and after."""

    original_data = wfi_model.data.copy()
    original_err = wfi_model.err.copy()

    step = SourceCatalogStep(
        snr_threshold=0.5,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
    )

    step.run(wfi_model)

    assert_allclose(original_data, wfi_model.data, atol=5.0e-7)
    assert_allclose(original_err, wfi_model.err, 5.0e-7)
