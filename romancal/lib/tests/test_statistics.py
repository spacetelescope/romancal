import numpy as np
import pytest
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel

from romancal.lib.statistics import populate_statistics


@pytest.mark.parametrize(
    "model_class",
    [ImageModel, MosaicModel],
)
@pytest.mark.parametrize(
    "data_val, dq_val, expected_frac",
    [
        # all pixels good (dq=0)
        (10.0, 0, 1.0),
        # all pixels flagged do_not_use (fraction should be 0.0)
        (10.0, pixel.DO_NOT_USE, 0.0),
        # half pixels flagged do_not_use (0.5)
        (5.0, "half_bad", 0.5),
        # other flags set (saturated), but not do_not_use (fraction should be 1.0)
        (10.0, pixel.SATURATED | pixel.HOT, 1.0),
    ],
)
def test_populate_statistics(model_class, data_val, dq_val, expected_frac):
    """Test statistics population."""
    shape = (10, 10)

    model = model_class.create_minimal()

    model.data = np.full(shape, data_val, dtype=np.float32)

    if dq_val == "half_bad":
        dq = np.zeros(shape, dtype=np.uint32)
        dq.view().flat[0:50] = pixel.DO_NOT_USE
        model.dq = dq
    else:
        model.dq = np.full(shape, dq_val, dtype=np.uint32)

    # inject a nan to verify nan-resistance (nanmedian and mad_std)
    model.data[0, 0] = np.nan

    populate_statistics(model)

    assert model.meta.statistics.image_median == data_val
    assert model.meta.statistics.image_rms == 0.0
    assert model.meta.statistics.good_pixel_fraction == pytest.approx(
        expected_frac, abs=1e-6
    )
    # check the placeholder logic (-1.0)
    assert model.meta.statistics.zodiacal_light == -1.0


def test_statistics_graceful_exit_no_data():
    """Ensure meta.statistics is created with defaults if data is None."""

    model = ImageModel.create_minimal()
    model.data = None

    populate_statistics(model)

    # 1. Check that meta.statistics WAS created (as per the first 'if' block)
    assert hasattr(model.meta, "statistics")

    # 2. Verify that default values were still assigned
    # Since model.data was None, the stats remain at their initialized defaults
    assert model.meta.statistics.zodiacal_light == -1.0
    assert np.isnan(model.meta.statistics.image_median)
    assert np.isnan(model.meta.statistics.image_rms)
    assert np.isnan(model.meta.statistics.good_pixel_fraction)
