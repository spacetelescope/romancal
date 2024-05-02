import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from numpy.testing import assert_allclose
from photutils.segmentation import SegmentationImage
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.maker_utils import mk_level2_image, mk_level3_mosaic

from romancal.source_catalog.reference_data import ReferenceData
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.source_catalog_step import SourceCatalogStep


def make_test_image():
    g1 = Gaussian2D(121.0, 11.1, 12.2, 1.5, 1.5)
    g2 = Gaussian2D(70, 65, 18, 9.2, 4.5)
    g3 = Gaussian2D(111.0, 41, 42.7, 8.0, 3.0, theta=30 * u.deg)
    g4 = Gaussian2D(81.0, 17, 52.7, 4, 2, theta=102 * u.deg)
    g5 = Gaussian2D(107.0, 65, 71, 12, 2, theta=142 * u.deg)
    g6 = Gaussian2D(50, 20, 80, 2.1, 2.1)
    g7 = Gaussian2D(97.0, 85, 87.3, 4, 2, theta=-30 * u.deg)

    yy, xx = np.mgrid[0:101, 0:101]
    data = (
        g1(xx, yy)
        + g2(xx, yy)
        + g3(xx, yy)
        + g4(xx, yy)
        + g5(xx, yy)
        + g6(xx, yy)
        + g7(xx, yy)
    ).value

    y0 = 2
    x0 = 90
    dd = 5
    value = -20
    y1 = y0 + dd
    x1 = x0 + dd - 1
    x2 = x0 + dd
    x3 = x0 + 2 * dd - 1
    data[y0:y1, x0] = value
    data[y0, x0:x1] = value
    data[y0:y1, x2] = value
    data[y0, x2:x3] = value
    data[y0 + dd // 2, x2:x3] = value
    data[y1 - 1, x2:x3] = value
    data[y0 + 1 : y1 - 1, x3] = value

    rng = np.random.default_rng(seed=123)
    noise = rng.normal(0, 2.5, size=data.shape)
    data += noise
    err = data / 10.0

    return data, err


@pytest.fixture
def mosaic_model():
    wfi_mosaic = mk_level3_mosaic()
    model = MosaicModel(wfi_mosaic)
    data, err = make_test_image()
    units = u.MJy / u.sr
    data <<= units
    err <<= units
    model.data = data
    model.err = err
    return model


@pytest.fixture
def image_model():
    wfi_image = mk_level2_image()
    model = ImageModel(wfi_image)
    data, err = make_test_image()
    units = u.DN / u.s
    data <<= units
    err <<= units
    model.data = data
    model.err = err
    model.meta.photometry.conversion_megajanskys = 0.3324 * u.MJy / u.sr
    return model


@pytest.mark.parametrize(
    "snr_threshold, npixels, nsources, save_results",
    (
        (3, 10, 7, True),
        (3, 50, 5, False),
        (10, 10, 7, False),
        (20, 10, 5, False),
        (25, 10, 3, False),
        (35, 10, 1, False),
        (50, 10, 0, False),
    ),
)
def test_l2_source_catalog(image_model, snr_threshold, npixels, nsources, save_results):
    step = SourceCatalogStep(
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )
    result = step.run(image_model)
    cat = result.source_catalog

    assert isinstance(cat, Table)

    columns = [
        "label",
        "xcentroid",
        "ycentroid",
        "ra_centroid",
        "dec_centroid",
        "aper_bkg_flux",
        "aper_bkg_flux_err",
        "aper30_flux",
        "aper30_flux_err",
        "aper50_flux",
        "aper50_flux_err",
        "aper70_flux",
        "aper70_flux_err",
        "aper_total_flux",
        "aper_total_flux_err",
        "CI_50_30",
        "CI_70_50",
        "CI_70_30",
        "is_extended",
        "sharpness",
        "roundness",
        "nn_label",
        "nn_dist",
        "isophotal_flux",
        "isophotal_flux_err",
        "isophotal_area",
        "semimajor_sigma",
        "semiminor_sigma",
        "ellipticity",
        "orientation",
        "sky_orientation",
        "ra_bbox_ll",
        "dec_bbox_ll",
        "ra_bbox_ul",
        "dec_bbox_ul",
        "ra_bbox_lr",
        "dec_bbox_lr",
        "ra_bbox_ur",
        "dec_bbox_ur",
    ]

    assert len(cat) == nsources

    if nsources > 0:
        for col in columns:
            assert col in cat.colnames
        assert np.min(cat["xcentroid"]) > 0.0
        assert np.min(cat["ycentroid"]) > 0.0
        assert np.max(cat["xcentroid"]) < 100.0
        assert np.max(cat["ycentroid"]) < 100.0


@pytest.mark.parametrize(
    "snr_threshold, npixels, nsources, save_results",
    (
        (3, 10, 7, True),
        (3, 50, 5, False),
        (10, 10, 7, False),
        (20, 10, 5, False),
        (25, 10, 3, False),
        (35, 10, 1, False),
        (50, 10, 0, False),
    ),
)
def test_l3_source_catalog(
    mosaic_model, snr_threshold, npixels, nsources, save_results
):
    step = SourceCatalogStep(
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )
    result = step.run(mosaic_model)
    cat = result.source_catalog

    assert isinstance(cat, Table)

    columns = [
        "label",
        "xcentroid",
        "ycentroid",
        "ra_centroid",
        "dec_centroid",
        "aper_bkg_flux",
        "aper_bkg_flux_err",
        "aper30_flux",
        "aper30_flux_err",
        "aper50_flux",
        "aper50_flux_err",
        "aper70_flux",
        "aper70_flux_err",
        "aper_total_flux",
        "aper_total_flux_err",
        "CI_50_30",
        "CI_70_50",
        "CI_70_30",
        "is_extended",
        "sharpness",
        "roundness",
        "nn_label",
        "nn_dist",
        "isophotal_flux",
        "isophotal_flux_err",
        "isophotal_area",
        "semimajor_sigma",
        "semiminor_sigma",
        "ellipticity",
        "orientation",
        "sky_orientation",
        "ra_bbox_ll",
        "dec_bbox_ll",
        "ra_bbox_ul",
        "dec_bbox_ul",
        "ra_bbox_lr",
        "dec_bbox_lr",
        "ra_bbox_ur",
        "dec_bbox_ur",
    ]

    assert len(cat) == nsources

    if nsources > 0:
        for col in columns:
            assert col in cat.colnames
        assert np.min(cat["xcentroid"]) > 0.0
        assert np.min(cat["ycentroid"]) > 0.0
        assert np.max(cat["xcentroid"]) < 100.0
        assert np.max(cat["ycentroid"]) < 100.0


def test_background(mosaic_model):
    """
    Test background fallback when Background2D fails.
    """
    step = SourceCatalogStep(
        bkg_boxsize=1000,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=25,
    )
    result = step.run(mosaic_model)
    cat = result.source_catalog

    assert isinstance(cat, Table)


def test_l2_input_model_unchanged(image_model):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    original_data = image_model.data.copy()
    original_err = image_model.err.copy()

    step = SourceCatalogStep(
        snr_threshold=0.5,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
    )
    step.run(image_model)

    assert_allclose(original_data, image_model.data, atol=5.0e-7)
    assert_allclose(original_err, image_model.err, atol=5.0e-7)


def test_l3_input_model_unchanged(mosaic_model):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    original_data = mosaic_model.data.copy()
    original_err = mosaic_model.err.copy()

    step = SourceCatalogStep(
        snr_threshold=0.5,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
    )
    step.run(mosaic_model)

    assert_allclose(original_data, mosaic_model.data, atol=5.0e-7)
    assert_allclose(original_err, mosaic_model.err, atol=5.0e-7)


def test_inputs(mosaic_model):
    with pytest.raises(ValueError):
        ReferenceData(np.ones((3, 3)), (30, 50, 70))
    with pytest.raises(ValueError):
        aperture_ee = (70, 50, 30)
        ReferenceData(mosaic_model, aperture_ee)
    with pytest.raises(ValueError):
        aperture_ee = (30, 50)
        ReferenceData(mosaic_model, aperture_ee)
    with pytest.raises(ValueError):
        aperture_ee = (-1, 50, 70)
        ReferenceData(mosaic_model, aperture_ee)
    with pytest.raises(ValueError):
        aperture_ee = (40, 70, 150)
        ReferenceData(mosaic_model, aperture_ee)

    data = np.ones((3, 3), dtype=int)
    data[1, 1] = 1
    segm = SegmentationImage(data)
    cdata = np.ones((3, 3))
    aper_params = {}
    ci_thresh = 100.0
    with pytest.raises(ValueError):
        RomanSourceCatalog(np.ones((3, 3)), segm, cdata, aper_params, ci_thresh, 2.0)

    with pytest.raises(ValueError):
        RomanSourceCatalog(mosaic_model, segm, cdata, aper_params, (1.0, 2.0, 3.0), 2.0)
