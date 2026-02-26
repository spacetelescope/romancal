from copy import deepcopy
from pathlib import Path
from re import match

import astropy.units as u
import numpy as np
import pyarrow
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from astropy.time import Time
from roman_datamodels import datamodels as rdm
from roman_datamodels.datamodels import MosaicModel, MultibandSegmentationMapModel

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog import MultibandCatalogStep
from romancal.multiband_catalog.multiband_catalog import match_recovered_sources
from romancal.skycell.tests.test_skycell_match import mk_gwcs


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
    ).value.astype("float32")

    rng = np.random.default_rng(seed=42)
    noise_scale = 2.5
    noise = rng.normal(0, noise_scale, size=data.shape)
    data += noise
    err = np.zeros_like(data) + noise_scale

    return data, err


@pytest.fixture
def mosaic_model(shape=(101, 101)):
    model = MosaicModel.create_fake_data(shape=shape)
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.var_rnoise = err**2
    model.weight = 1.0 / err
    model.meta.instrument.optical_element = "F184"
    model.meta.coadd_info.time_first = Time("2027-01-01T00:00:00")
    model.meta.wcsinfo.pixel_scale = 0.11 / 3600  # degrees

    model.meta.wcsinfo.ra_ref = 270.0  # degrees
    model.meta.wcsinfo.dec_ref = 66.0  # degrees
    model.meta.wcsinfo.roll_ref = 0.0  # degrees
    model.meta.coadd_info.exposure_time = 300  # seconds

    model.meta.resample.pixfrac = 0.5
    model.meta.data_release_id = "r1"
    return model


@pytest.fixture
def library_model(mosaic_model):
    model2 = mosaic_model.copy()
    model2.meta.instrument.optical_element = "F158"
    return ModelLibrary([mosaic_model, model2])


@pytest.fixture
def library_model_all_nan(mosaic_model):
    model1 = mosaic_model.copy()
    model1.data[:] = np.nan
    model2 = mosaic_model.copy()
    model2.data[:] = np.nan
    model2.meta.instrument.optical_element = "F158"
    return ModelLibrary([model1, model2])


def shared_tests(
    result, cat, library_model, save_results, function_jail, shape=(101, 101)
):
    with library_model:
        input_model = library_model.borrow(0)
        assert result.meta.data_release_id == input_model.meta.data_release_id
        library_model.shelve(input_model, modify=False)

    assert len(cat.meta["aperture_radii"]["circle_pix"]) > 0
    assert sum(
        1 for name in cat.colnames if match(r"^aper\d+_f158_flux$", name)
    ) == len(cat.meta["aperture_radii"]["circle_pix"])
    assert sum(
        1 for name in cat.colnames if match(r"^aper\d+_f184_flux$", name)
    ) == len(cat.meta["aperture_radii"]["circle_pix"])
    assert "ee_fractions" in cat.meta
    assert isinstance(cat.meta["ee_fractions"], dict)
    assert len(cat.meta["ee_fractions"]) == 2
    assert "f158" in cat.meta["ee_fractions"]
    assert "f184" in cat.meta["ee_fractions"]
    for value in cat.meta["ee_fractions"].values():
        assert len(value) == len(cat.meta["aperture_radii"]["circle_pix"])

    if len(cat) > 0:
        assert np.min(cat["x_centroid"]) > 0.0
        assert np.min(cat["y_centroid"]) > 0.0
        assert np.max(cat["x_centroid"]) < shape[0]
        assert np.max(cat["y_centroid"]) < shape[1]

        for colname in cat.colnames:
            if (
                "flux" in colname
                and "fluxfrac" not in colname
                and "aper_bkg" not in colname
            ):
                assert cat[colname].unit == "nJy"
                assert "f158" in colname or "f184" in colname
            if colname.endswith("_flux"):
                assert f"{colname}_err" in cat.colnames

    if save_results:
        filepath = Path(function_jail / f"{result.meta.filename}_cat.parquet")
        assert filepath.exists()
        tbl = pyarrow.parquet.read_table(filepath)
        assert isinstance(tbl, pyarrow.Table)

        filepath = Path(function_jail / f"{result.meta.filename}_segm.asdf")
        assert filepath.exists()
        assert isinstance(rdm.open(filepath), MultibandSegmentationMapModel)


@pytest.mark.parametrize("fit_psf", (True, False))
@pytest.mark.parametrize(
    "snr_threshold, npixels, save_results",
    (
        (3, 10, True),
        (7, 10, False),
    ),
)
def test_multiband_catalog(
    library_model, fit_psf, snr_threshold, npixels, save_results, function_jail
):
    step = MultibandCatalogStep()

    result = step.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=snr_threshold,
        npixels=npixels,
        fit_psf=fit_psf,
        save_results=save_results,
        deblend=True,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7

    shared_tests(result, cat, library_model, save_results, function_jail)


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_no_detections(library_model, save_results, function_jail):
    step = MultibandCatalogStep()

    result = step.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=1000,  # high threshold to ensure no detections
        npixels=10,
        fit_psf=False,
        save_results=save_results,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 0


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_invalid_inputs(
    library_model_all_nan, save_results, function_jail
):
    step = MultibandCatalogStep()

    result = step.call(
        library_model_all_nan,
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=False,
        save_results=save_results,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 0


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_some_invalid_inputs(
    library_model, save_results, function_jail
):
    # Modify the first model in the library to have all NaN values
    with library_model:
        model = library_model.borrow(0)  # f184 model
        model.data[:] = np.nan
        model.err[:] = np.nan
        model.var_rnoise[:] = np.nan
        library_model.shelve(model, modify=True)

    step = MultibandCatalogStep()

    result = step.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=False,
        save_results=save_results,
        deblend=True,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7
    assert "segment_f158_flux" in cat.colnames
    assert "segment_f184_flux" in cat.colnames
    assert np.all(np.isnan(cat["segment_f184_flux"]))
    assert np.all(np.isnan(cat["segment_f184_flux_err"]))


def make_si_test_image():
    g1 = Gaussian2D(60.5, 11, 12, 1.5, 1.5)
    g2 = Gaussian2D(35, 65, 18, 9.2, 4.5)
    g3 = Gaussian2D(55.5, 41, 43, 8.0, 3.0, theta=30 * u.deg)
    g4 = Gaussian2D(40.5, 17, 53, 4, 2, theta=102 * u.deg)
    g5 = Gaussian2D(53.5, 65, 71, 12, 2, theta=142 * u.deg)
    g6 = Gaussian2D(25, 20, 80, 2.1, 2.1)
    g7 = Gaussian2D(48.5, 85, 88, 4, 2, theta=-30 * u.deg)

    yy, xx = np.mgrid[0:100, 0:100]
    smalldata = np.zeros(shape=(len(yy), len(xx)))
    smalldata = (
        g1(xx, yy)
        + g2(xx, yy)
        + g3(xx, yy)
        + g4(xx, yy)
        + g5(xx, yy)
        + g6(xx, yy)
        + g7(xx, yy)
    ).value.astype("float32")

    data = np.zeros(shape=(500, 500))
    data = np.tile(smalldata, (5, 5))

    rng = np.random.default_rng(seed=42)
    noise_scale = 0.01 * 0.2
    noise = rng.normal(0, noise_scale, size=data.shape)
    data += noise
    err = np.zeros_like(data) + noise_scale

    return data, err


@pytest.fixture
def mosaic_si_model(shape=(500, 500)):
    model = MosaicModel.create_fake_data(shape=shape)
    data, err = make_si_test_image()
    model.data = data
    model.err = err
    model.var_rnoise = err**2
    model.var_poisson = err**2
    model.weight = 1.0 / err
    model.meta.instrument.optical_element = "F184"
    model.meta.coadd_info.time_first = Time("2027-01-01T00:00:00")
    model.meta.wcsinfo.pixel_scale = 0.11 / 3600  # degrees

    model.meta.wcsinfo.ra_ref = 270.0  # degrees
    model.meta.wcsinfo.dec_ref = 66.0  # degrees
    model.meta.wcsinfo.roll_ref = 0.0  # degrees
    model.meta.coadd_info.exposure_time = 1  # seconds

    model.meta.resample.pixfrac = 0.5
    model.meta.data_release_id = "r1"

    # Create WCS
    model.meta.wcs = mk_gwcs(
        model.meta.wcsinfo.ra_ref,
        model.meta.wcsinfo.dec_ref,
        model.meta.wcsinfo.roll_ref,
        bounding_box=((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5)),
        shape=shape,
    )

    return model


@pytest.fixture
def library_model2(mosaic_si_model):
    si_model2 = deepcopy(mosaic_si_model)
    si_model2.meta.instrument.optical_element = "F158"
    return ModelLibrary([mosaic_si_model, si_model2])


@pytest.mark.parametrize("fit_psf", (True, False))
@pytest.mark.parametrize(
    "snr_threshold, npixels, save_results",
    (
        (3, 10, False),
        (2, 4, True),
    ),
)
def test_multiband_source_injection_catalog(
    library_model2, fit_psf, snr_threshold, npixels, save_results, function_jail
):
    step = MultibandCatalogStep()

    result = step.call(
        library_model2,
        bkg_boxsize=50,
        snr_threshold=snr_threshold,
        npixels=npixels,
        fit_psf=fit_psf,
        deblend=True,
        inject_sources=True,
        inject_seed=50,
        save_results=save_results,
        save_debug_info=True,
    )

    # Original objects
    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 175

    # Ensure all original objects found in the proper location
    cat["x_mod"] = np.mod(np.round(cat["x_centroid"]), 100)
    cat["y_mod"] = np.mod(np.round(cat["y_centroid"]), 100)
    gocj_locs = [
        (11, 12),
        (65, 18),
        (41, 43),
        (17, 53),
        (65, 71),
        (20, 80),
        (85, 88),
    ]
    for modx, mody in cat[["x_mod", "y_mod"]]:
        assert (modx, mody) in gocj_locs

    # Source injected and original images
    si_cat = result.source_injection_catalog
    assert isinstance(si_cat, Table)
    assert len(si_cat) >= len(cat)

    if save_results:
        filepath = Path(function_jail / f"{result.meta.filename}_segm.asdf")
        assert filepath.exists()
        segm_mod = rdm.open(filepath)
        assert isinstance(segm_mod, MultibandSegmentationMapModel)
        assert segm_mod.data.shape == segm_mod.si_data.shape
        assert isinstance(segm_mod.injected_sources, Table)
        assert len(segm_mod.injected_sources[0]) <= len(si_cat)

        assert np.count_nonzero(
            segm_mod.recovered_sources["best_injected_index"] != -1
        ) >= (400 / 2)

    # Old lines from other MBC tests
    shared_tests(
        result,
        cat,
        library_model2,
        test_multiband_catalog,
        function_jail,
        shape=(5000, 5000),
    )


def test_match_recovered_sources():
    # Create minimal source catalog data
    shape = (2000, 2000)
    orig_xy = [[510, 510], [400, 400], [1000, 1000], [450, 1450], [100, 100]]
    injected_xy = [[500, 500], [500, 1500], [1500, 500], [1500, 1500]]
    injected_hlr = [2.0, 12.0, 5.0, 1.0]
    si_cat_xy = [
        [510.1, 510.1],
        [400.1, 400.1],
        [1000.1, 1000.1],
        [450.1, 1450.1],
        [100.1, 100.1],
        [500.1, 500.1],
        [500.1, 1500.1],
        [1500.1, 500.1],
        [1500.1, 1500.1],
    ]

    # Make wcs object
    wcsobj = mk_gwcs(
        270.0,  # degrees
        66.0,  # degrees
        0.0,  # degrees
        bounding_box=((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5)),
        shape=shape,
    )

    # Create source catalogs
    orig_table = Table()
    injected_table = Table()
    si_cat_table = Table()

    # Conert x, y to ra & dec
    orig_table["x"], orig_table["y"] = zip(*orig_xy, strict=True)
    orig_table["ra"], orig_table["dec"] = wcsobj.pixel_to_world_values(
        np.array(orig_table["x"]), np.array(orig_table["y"])
    )

    injected_table["x"], injected_table["y"] = zip(*injected_xy, strict=True)
    injected_table["ra"], injected_table["dec"] = wcsobj.pixel_to_world_values(
        np.array(injected_table["x"]), np.array(injected_table["y"])
    )

    si_cat_table["x"], si_cat_table["y"] = zip(*si_cat_xy, strict=True)
    si_cat_table["ra"], si_cat_table["dec"] = wcsobj.pixel_to_world_values(
        np.array(si_cat_table["x"]), np.array(si_cat_table["y"])
    )

    # Additional table columns
    injected_table["half_light_radius"] = injected_hlr
    injected_table["half_light_radius"].unit = "arcsec"
    orig_table["empty"] = 0
    si_cat_table["one"] = 1

    # Conduct matching
    rec_table = match_recovered_sources(orig_table, injected_table, si_cat_table)

    # Test that distant sources were dropped
    assert len(rec_table) < len(si_cat_table)

    # Test that grid sources were found
    assert np.count_nonzero(rec_table["best_injected_index"] != -1) == len(
        injected_table
    )

    # Test columns included or excluded as expected
    assert "one" in rec_table.colnames
    assert "empty" not in rec_table.colnames
