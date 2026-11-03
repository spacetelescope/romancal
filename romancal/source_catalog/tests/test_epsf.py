"""Unit tests for the ePSF star selection."""

import numpy as np
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.dqflags import pixel

from romancal.source_catalog import _epsf

SHAPE = (200, 200)


@pytest.fixture
def image_and_catalog():
    """Isolated stars on a grid plus a few special cases.

    Labels 1-25 are isolated, unflagged stars. Label 26 is saturated,
    27 has a DO_NOT_USE core, 28 has a DO_NOT_USE wing pixel, 29 is
    faint, and 30/31 are a blended pair.
    """
    rng = np.random.default_rng(42)
    yy, xx = np.mgrid[: SHAPE[0], : SHAPE[1]]
    rows = []
    for i, (y, x) in enumerate(
        (y, x) for y in range(20, 181, 40) for x in range(20, 181, 40)
    ):
        rows.append((i + 1, x + rng.uniform(-0.5, 0.5), y + rng.uniform(-0.5, 0.5), 1e4))
    rows += [
        (26, 40.3, 60.2, 1e4),
        (27, 80.1, 60.4, 1e4),
        (28, 120.4, 60.1, 1e4),
        (29, 160.2, 60.3, 30.0),
        (30, 40.2, 140.1, 1e4),
        (31, 43.2, 140.3, 5e3),
    ]
    cat = Table(rows=rows, names=["label", "x", "y", "flux"])

    sigma = 1.0
    data = np.zeros(SHAPE)
    for row in cat:
        amp = row["flux"] / (2 * np.pi * sigma**2)
        data += Gaussian2D(amp, row["x"], row["y"], sigma, sigma)(xx, yy)
    err = np.full(SHAPE, 1.0)
    data += rng.normal(0, 1.0, SHAPE)

    dq = np.zeros(SHAPE, dtype=np.uint32)
    dq[60, 40] = pixel.SATURATED
    dq[60, 80] = pixel.DO_NOT_USE
    dq[63, 120] = pixel.DO_NOT_USE
    data[63, 120] = np.nan

    model = ImageModel.create_fake_data(shape=SHAPE)
    model.data = data.astype(np.float32)
    model.err = err.astype(np.float32)
    model.var_poisson = np.zeros(SHAPE, dtype=np.float32)
    model.var_rnoise = (err**2).astype(np.float32)
    model.dq = dq

    snr = cat["flux"] / (np.sqrt(4 * np.pi) * sigma)
    catalog = Table(
        {
            "label": cat["label"],
            "x_psf": cat["x"],
            "y_psf": cat["y"],
            "x_centroid": cat["x"],
            "y_centroid": cat["y"],
            "x_psf_err": np.full(len(cat), 0.01),
            "y_psf_err": np.full(len(cat), 0.01),
            "psf_flux": cat["flux"],
            "psf_flux_err": cat["flux"] / snr,
            "psf_flags": np.zeros(len(cat), dtype=int),
            "is_extended": np.zeros(len(cat), dtype=bool),
            "segment_flux": cat["flux"],
        }
    )
    return model, catalog


def picked_labels(model, catalog, **kwargs):
    stars, _ = _epsf.pick_stars(model, catalog, max_psferr=1.0, **kwargs)
    return {int(star.id_label) for star in stars}


def test_pick_stars_default(image_and_catalog):
    """Isolated stars are kept; special cases are handled."""
    labels = picked_labels(*image_and_catalog)
    assert set(range(1, 26)) <= labels
    assert 26 not in labels  # saturated
    assert 27 not in labels  # DO_NOT_USE core
    assert 28 in labels  # DO_NOT_USE wing pixel is masked, not rejected
    assert 29 not in labels  # faint
    assert not {30, 31} & labels  # blended


def test_pick_stars_masks_wing_pixel(image_and_catalog):
    """A bad wing pixel is masked in the stamp."""
    stars, _ = _epsf.pick_stars(*image_and_catalog, max_psferr=1.0)
    (star,) = [s for s in stars if s.id_label == 28]
    assert star.mask.sum() == 1


@pytest.mark.parametrize(
    "kwargs, label",
    [
        ({"reject_saturated": False}, 26),
        ({"reject_core_dq": None}, 27),
        ({"min_snr": None}, 29),
        ({"max_contamination": None}, 30),
    ],
)
def test_pick_stars_options(image_and_catalog, kwargs, label):
    """Each rejection can be turned off."""
    assert label in picked_labels(*image_and_catalog, **kwargs)


@pytest.mark.parametrize("flag, kept", [(1, True), (8, False)])
def test_pick_stars_psf_flags(image_and_catalog, flag, kept):
    """Only the selected psf_flags bits reject a star."""
    model, catalog = image_and_catalog
    catalog["psf_flags"][0] = flag
    assert (1 in picked_labels(model, catalog)) == kept


def test_background_noise(image_and_catalog):
    """Noise is read noise plus sky, with a fallback without var_rnoise."""
    model, _ = image_and_catalog
    model.var_poisson[:] = 0.5
    model.var_poisson[100, 100] = 100.0  # bright pixel
    model.var_rnoise[:] = 2.0
    model.err = np.sqrt(model.var_poisson + model.var_rnoise)
    np.testing.assert_allclose(_epsf.background_noise(model), np.sqrt(2.5))
    delattr(model, "var_rnoise")
    np.testing.assert_allclose(_epsf.background_noise(model), np.sqrt(2.5))
