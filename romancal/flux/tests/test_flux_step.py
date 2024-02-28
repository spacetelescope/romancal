"""Unit-like tests related to FluxStep"""
import numpy as np
from numpy.random import poisson
import pytest

import astropy.units as u

from roman_datamodels import datamodels, maker_utils
from romancal.flux import FluxStep


@pytest.mark.parametrize(
    'attr, factor',
    [
        ('data', 1),
        ('var_rnoise', 2),
        ('var_poisson', 2),
        ('var_flat', 2),
    ]
)
def test_attributes(flux_step, attr, factor):
    """Test that the attribute has been scaled by the right factor"""
    original, result = flux_step

    c_mj = original.meta.photometry.conversion_megajanskys
    scale = c_mj**factor
    original_value = getattr(original, attr)
    result_value = getattr(result, attr)

    assert np.allclose(original_value * scale, result_value)


# ########
# Fixtures
# ########
@pytest.fixture(scope='module')
def flux_step(input):
    """Execute FluxStep on given input

    Parameters
    ----------
    input : str, `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.container.ModelContainer`

    Returns
    -------
    original, result : DataModel or ModelContainer, DataModel or ModelContainer
    """

    # Copy input because flux operates in-place
    original = input.copy()

    # Perform step
    result = FluxStep.call(input)

    # That's all folks
    return original, result


@pytest.fixture(scope='module')
def image_model():
    """Product a basic ImageModel"""
    # Create a random image and specify a conversion.
    rng = np.random.default_rng()
    shape = (10, 10)
    image_model = maker_utils.mk_datamodel(datamodels.ImageModel, shape=shape)
    image_model.data = u.Quantity(
        rng.poisson(2.5, size=shape).astype(np.float32), u.electron / u.s, dtype=np.float32)
    image_model.var_rnoise = u.Quantity(
        rng.normal(1, 0.05, size=shape).astype(np.float32), u.electron**2 / u.s**2, dtype=np.float32)
    image_model.var_poisson = u.Quantity(
        rng.poisson(1, size=shape).astype(np.float32), u.electron**2 / u.s**2, dtype=np.float32)
    image_model.var_flat = u.Quantity(
        rng.uniform(0, 1, size=shape).astype(np.float32), u.electron**2 / u.s**2, dtype=np.float32)
    image_model.meta.photometry.conversion_megajanskys = 2.0 * u.MJy / u.sr

    return image_model

@pytest.fixture(scope='module')
def input(image_model):
    """Setup inputs to the FluxStep"""

    # First just setup the basic model
    yield image_model
