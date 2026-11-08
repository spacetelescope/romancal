"""Tests for resample utilities."""

import numpy as np
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.dqflags import pixel

from romancal.resample._utils import compute_var_sky

SHAPE = (20, 20)


def make_model(data):
    """Model with inverse gain 0.5 and read noise variance 3."""
    model = ImageModel.create_fake_data(shape=SHAPE)
    model.data = data.astype(np.float32)
    model.dq = np.zeros(SHAPE, dtype=np.uint32)
    model.var_poisson = (0.5 * np.clip(data, 0, None)).astype(np.float32)
    model.err = np.sqrt(model.var_poisson + 3.0).astype(np.float32)
    for name in ("var_rnoise", "var_flat", "var_dark"):
        if hasattr(model, name):
            delattr(model, name)
    return model


def test_compute_var_sky():
    """Sky variance is independent of flux, including at low and NaN data."""
    data = np.full(SHAPE, 2.0)
    data[0, 0] = 1000.0  # bright
    data[0, 1] = -1.0  # negative; var_poisson clipped at zero
    data[0, 2] = 0.01  # near zero
    data[0, 3] = np.nan  # without DO_NOT_USE
    var_sky = compute_var_sky(make_model(data))
    expected = np.full(SHAPE, 3.0 + 0.5 * 2.0)
    expected[0, 3] = np.nan
    np.testing.assert_allclose(var_sky, expected, rtol=1e-3)


def test_compute_var_sky_no_sky():
    """Without a positive sky level, the sky variance is the read noise."""
    data = np.zeros(SHAPE)
    model = make_model(data)
    model.dq[:] = 0
    model.dq[0, 0] = pixel.DO_NOT_USE
    np.testing.assert_allclose(compute_var_sky(model), 3.0, rtol=1e-3)
