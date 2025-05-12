import os
from pathlib import Path

import astropy.units as u
import numpy as np
import pyarrow
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from roman_datamodels import datamodels as rdm
from roman_datamodels.datamodels import MosaicModel, MosaicSegmentationMapModel
from roman_datamodels.maker_utils import mk_level3_mosaic

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
    wfi_mosaic = mk_level3_mosaic(shape=(101, 101))
    model = MosaicModel(wfi_mosaic)
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.var_rnoise = err**2
    model.weight = 1.0 / err
    return model


@pytest.fixture
def library_model(mosaic_model):
    model2 = mosaic_model.copy()
    model2.meta.basic.optical_element = "F184"
    return ModelLibrary([mosaic_model, model2])


@pytest.mark.stpsf
@pytest.mark.parametrize(
    "snr_threshold, npixels, save_results",
    (
        (3, 10, True),
        (7, 10, False),
    ),
)
def test_multiband_catalog(
    library_model, snr_threshold, npixels, save_results, tmp_path
):
    os.chdir(tmp_path)
    step = MultibandCatalogStep()

    result = step.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=snr_threshold,
        npixels=npixels,
        fit_psf=False,
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
