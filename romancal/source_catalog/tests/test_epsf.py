"""Unit tests for the ePSF star selection."""

import numpy as np
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.dqflags import pixel

from romancal.source_catalog import _epsf
from romancal.source_catalog.psf import add_jitter

SHAPE = (200, 200)


@pytest.fixture
def image_and_catalog():
    """Isolated stars on a grid plus a few special cases.

    Labels 1-25 are isolated, unflagged stars. Label 26 is saturated,
    27 has a DO_NOT_USE core, 28 has a DO_NOT_USE wing pixel, 29 is
    faint, 30/31 are a blended pair, and 32 is a single hot pixel.
    """
    rng = np.random.default_rng(42)
    yy, xx = np.mgrid[: SHAPE[0], : SHAPE[1]]
    rows = []  # label, x, y, flux
    for i, (y, x) in enumerate(
        (y, x) for y in range(20, 181, 40) for x in range(20, 181, 40)
    ):  # grid of 5x5 bright stars in a 200x200 image
        rows.append(
            (i + 1, x + rng.uniform(-0.5, 0.5), y + rng.uniform(-0.5, 0.5), 1e4)
        )
    rows += [
        (26, 40.3, 60.2, 1e4),
        (27, 80.1, 60.4, 1e4),
        (28, 120.4, 60.1, 1e4),
        (29, 160.2, 60.3, 30.0),
        (30, 40.2, 140.1, 1e4),
        (31, 43.2, 140.3, 5e3),
        (32, 120.0, 140.0, 1e4),
    ]  # extra problem cases
    cat = Table(rows=rows, names=["label", "x", "y", "flux"])

    sigma = 1.0
    data = np.zeros(SHAPE)
    for row in cat:
        amp = row["flux"] / (2 * np.pi * sigma**2)
        # simple Gaussian PSF profile
        data += Gaussian2D(amp, row["x"], row["y"], sigma, sigma)(xx, yy)
    err = np.full(SHAPE, 1.0)
    data += rng.normal(0, 1.0, SHAPE)

    # Label 32 is a single hot pixel, not a star.
    data -= Gaussian2D(
        cat["flux"][-1] / (2 * np.pi * sigma**2), 120.0, 140.0, sigma, sigma
    )(xx, yy)
    data[140, 120] += 1e4

    # add some saturated / DO_NOT_USE stars
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
            "psf_flux": cat["flux"],
            "psf_flux_err": cat["flux"] / snr,
            "psf_flags": np.zeros(len(cat), dtype=int),
            "is_extended": np.zeros(len(cat), dtype=bool),
            "kron_flux": cat["flux"],
            # The blended pair sees each other's flux without masking.
            "kron_nomask_flux": np.where(
                np.isin(cat["label"], [30, 31]), 2 * cat["flux"], cat["flux"]
            ),
            # The hot pixel (last row) is badly fit by a PSF.
            "psf_cfit": np.where(cat["label"] == 32, 2.0, 0.0),
        }
    )
    return model, catalog


# Test values for cut parameters
CUTS = {
    "stamp_size": 21,
    "reject_core_dq": pixel.DO_NOT_USE,
    "reject_saturated": True,
    "max_masked_fraction": 0.2,
    "max_psf_cfit": 0.2,
    "min_snr": 50.0,
    "max_kron_contamination": 0.5,
    "reject_extended": True,
    "reject_psf_flags": 8,
}


def picked_labels(model, catalog, **kwargs):
    stars = _epsf.pick_stars(model, catalog, **{**CUTS, **kwargs})
    return {int(star.id_label) for star in stars}


def test_pick_stars(image_and_catalog):
    """Isolated stars are kept; special cases are rejected."""
    labels = picked_labels(*image_and_catalog)
    assert labels == set(range(1, 26)) | {28}


def test_pick_stars_masks_wing_pixel(image_and_catalog):
    """A bad wing pixel is masked in the stamp, not rejected."""
    stars = _epsf.pick_stars(*image_and_catalog, **CUTS)
    star = next(s for s in stars if s.id_label == 28)
    assert star.mask.sum() == 1


@pytest.mark.parametrize(
    "kwargs, label",
    [  # test that turning off various cuts changes whether appropriate stars are rejected
        ({"reject_saturated": False}, 26),
        ({"reject_core_dq": None}, 27),
        ({"min_snr": None}, 29),
        ({"max_kron_contamination": None}, 30),
        ({"max_psf_cfit": None}, 32),
        ({"reject_extended": False}, 1),
        ({"reject_psf_flags": None}, 2),
    ],
)
def test_pick_stars_options(image_and_catalog, kwargs, label):
    """Each cut can be turned off."""
    model, catalog = image_and_catalog
    catalog["is_extended"][0] = True
    catalog["psf_flags"][1] = 8
    assert label not in picked_labels(model, catalog)
    assert label in picked_labels(model, catalog, **kwargs)


@pytest.mark.parametrize("flag, kept", [(1, True), (8, False), (9, False)])
def test_pick_stars_psf_flags(image_and_catalog, flag, kept):
    """Only the selected psf_flags bits reject a star."""
    model, catalog = image_and_catalog
    catalog["psf_flags"][0] = flag
    assert (1 in picked_labels(model, catalog, reject_psf_flags=8)) == kept


def test_epsf_from_model_jitter(image_and_catalog):
    """The reference records the image jitter, so add_jitter is a no-op."""
    model, catalog = image_and_catalog
    model.meta.guide_star.jitter_major = 12.0
    model.meta.guide_star.jitter_minor = 6.0
    model.meta.guide_star.jitter_position_angle = 30.0
    epsf_model, _, _ = _epsf.epsf_from_model(model, catalog, kwargs_pick_stars=CUTS)
    assert epsf_model.meta.jitter_major == 12.0
    assert epsf_model.meta.jitter_minor == 6.0
    assert epsf_model.meta.jitter_position_angle == 30.0
    np.testing.assert_allclose(add_jitter(epsf_model, model), epsf_model.psf, atol=1e-6)
