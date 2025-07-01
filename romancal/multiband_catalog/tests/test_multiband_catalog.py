import os
from pathlib import Path

import astropy.units as u
import numpy as np
import pyarrow
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from astropy.time import Time
from roman_datamodels import datamodels as rdm
from roman_datamodels.datamodels import MosaicModel, MosaicSegmentationMapModel

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog import MultibandCatalogStep


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

    rng = np.random.default_rng(seed=123)
    noise_scale = 2.5
    noise = rng.normal(0, noise_scale, size=data.shape)
    data += noise
    err = np.zeros_like(data) + noise_scale

    return data, err


@pytest.fixture
def mosaic_model():
    model = MosaicModel.create_fake_data(shape=(101, 101))
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.var_rnoise = err**2
    model.weight = 1.0 / err
    model.meta.basic.optical_element = "F184"
    model.meta.basic.time_first_mjd = Time("2027-01-01T00:00:00").mjd
    model.meta.wcsinfo.pixel_scale = 0.11 / 3600  # degrees
    model.meta.resample.pixfrac = 0.5
    return model


@pytest.fixture
def library_model(mosaic_model):
    model2 = mosaic_model.copy()
    model2.meta.basic.optical_element = "F158"
    return ModelLibrary([mosaic_model, model2])


@pytest.fixture
def library_model_all_nan(mosaic_model):
    model1 = mosaic_model.copy()
    model1.data[:] = np.nan
    model2 = mosaic_model.copy()
    model2.data[:] = np.nan
    model2.meta.basic.optical_element = "F158"
    return ModelLibrary([model1, model2])


@pytest.mark.parametrize("fit_psf", (True, False))
@pytest.mark.parametrize(
    "snr_threshold, npixels, save_results",
    (
        (3, 10, True),
        (7, 10, False),
    ),
)
def test_multiband_catalog(
    library_model, fit_psf, snr_threshold, npixels, save_results, tmp_path
):
    os.chdir(tmp_path)
    step = MultibandCatalogStep()

    result = step.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=snr_threshold,
        npixels=npixels,
        fit_psf=fit_psf,
        save_results=save_results,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7

    if len(cat) > 0:
        assert np.min(cat["x_centroid"]) > 0.0
        assert np.min(cat["y_centroid"]) > 0.0
        assert np.max(cat["x_centroid"]) < 100.0
        assert np.max(cat["y_centroid"]) < 100.0

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
        filepath = Path(tmp_path / f"{result.meta.filename}_cat.parquet")
        assert filepath.exists()
        tbl = pyarrow.parquet.read_table(filepath)
        assert isinstance(tbl, pyarrow.Table)

        filepath = Path(tmp_path / f"{result.meta.filename}_segm.asdf")
        assert filepath.exists()
        assert isinstance(rdm.open(filepath), MosaicSegmentationMapModel)


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_no_detections(library_model, save_results, tmp_path):
    os.chdir(tmp_path)
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
    library_model_all_nan, save_results, tmp_path
):
    os.chdir(tmp_path)
    step = MultibandCatalogStep()

    result = step.call(
        library_model_all_nan,
        bkg_boxsize=50,
        snr_threshold=1000,  # high threshold to ensure no detections
        npixels=10,
        fit_psf=False,
        save_results=save_results,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 0
