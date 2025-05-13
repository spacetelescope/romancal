"""Unit-like tests related to FluxStep"""

import astropy.units as u
import numpy as np
import pytest
from roman_datamodels import datamodels, stnode

from romancal.datamodels import ModelLibrary
from romancal.flux import FluxStep


@pytest.mark.parametrize(
    "attr, factor",
    [
        ("data", 1),
        ("err", 1),
        ("var_rnoise", 2),
        ("var_poisson", 2),
        ("var_flat", 2),
    ],
)
def test_attributes(flux_step, attr, factor):
    """Test that the attribute has been scaled by the right factor"""
    original, result = flux_step
    c_unit = 1.0

    # Handle difference between just a single image and a list.
    if isinstance(original, datamodels.ImageModel):
        original_library = ModelLibrary([original])
        result_library = ModelLibrary([result])
    else:
        original_library = original
        result_library = result

    assert len(original_library) == len(result_library)
    with original_library, result_library:
        for i in range(len(original_library)):
            original_model = original_library.borrow(i)
            result_model = result_library.borrow(i)

            c_mj = original_model.meta.photometry.conversion_megajanskys
            scale = (c_mj * c_unit) ** factor
            original_value = getattr(original_model, attr)
            result_value = getattr(result_model, attr)

            assert np.allclose(original_value * scale, result_value)

            original_library.shelve(original_model, i, modify=False)
            result_library.shelve(result_model, i, modify=False)


# ########
# Fixtures
# ########
@pytest.fixture(scope="module", params=["input_imagemodel", "input_modellibrary"])
def flux_step(request):
    """Execute FluxStep on given input

    Parameters
    ----------
    input : str, `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.library.ModelLibrary`

    Returns
    -------
    original, result : DataModel or ModelLibrary, DataModel or ModelLibrary
    """
    init = request.getfixturevalue(request.param)
    if isinstance(init, ModelLibrary):
        models = []
        with init:
            for m in init:
                models.append(m.copy())
                init.shelve(m, modify=False)
        original = ModelLibrary(models)
    else:
        original = init.copy()

    # Perform step
    result = FluxStep.call(init)

    # That's all folks
    return original, result


@pytest.fixture(scope="module")
def image_model():
    """Product a basic ImageModel"""
    # Create a random image and specify a conversion
    rng = np.random.default_rng()
    shape = (10, 10)
    image_model = datamodels.ImageModel.fake_data(shape=shape)
    image_model.data = rng.poisson(2.5, size=shape).astype(np.float32)
    image_model.var_rnoise = rng.normal(1, 0.05, size=shape).astype(np.float32)
    image_model.var_poisson = rng.poisson(1, size=shape).astype(np.float32)
    image_model.var_flat = rng.uniform(0, 1, size=shape).astype(np.float32)
    image_model.meta.photometry.conversion_megajanskys = (2.0 * u.MJy / u.sr).value
    image_model.meta.cal_step = stnode.L2CalStep.fake_data()
    image_model.meta.cal_logs = stnode.CalLogs.fake_data()

    return image_model


@pytest.fixture(scope="module")
def input_imagemodel(image_model):
    """Provide a single ImageModel"""

    # First just setup the basic model
    return image_model.copy()


@pytest.fixture(scope="module")
def input_modellibrary(image_model):
    """Provide a ModelLibrary"""
    # Create and return a ModelLibrary
    image_model1 = image_model.copy()
    image_model2 = image_model.copy()
    image_model2.meta.photometry.conversion_megajanskys = (0.5 * u.MJy / u.sr).value
    container = ModelLibrary([image_model1, image_model2])
    return container
