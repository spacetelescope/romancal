import numpy as np
import pytest
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel

from romancal.lib.statistics import populate_statistics

SHAPE = (10, 10)


def _create_model(model_class, *, data, dq=None, err=None):
    model = model_class.create_minimal()
    model.data = np.asarray(data, dtype=np.float32)

    if dq is not None:
        model.dq = np.asarray(dq, dtype=np.uint32)

    model.err = None if err is None else np.asarray(err, dtype=np.float32)

    return model


def _assert_float_or_nan(actual, expected):
    if np.isnan(expected):
        assert np.isnan(actual)
    else:
        assert actual == pytest.approx(expected, abs=1e-6)


def _assert_model_specific_placeholders(model):
    if isinstance(model, MosaicModel):
        assert not hasattr(model.meta.statistics, "zodiacal_light")
    else:
        assert model.meta.statistics.zodiacal_light == -1.0


def _case_all_pixels_good():
    data = np.full(SHAPE, 10.0, dtype=np.float32)
    dq = np.zeros(SHAPE, dtype=np.uint32)
    return data, dq, None, 10.0, 0.0, 1.0


def _case_all_pixels_do_not_use():
    data = np.full(SHAPE, 10.0, dtype=np.float32)
    dq = np.full(SHAPE, pixel.DO_NOT_USE, dtype=np.uint32)
    return data, dq, None, np.nan, np.nan, 0.0


def _case_half_pixels_do_not_use():
    data = np.full(SHAPE, 5.0, dtype=np.float32)
    dq = np.zeros(SHAPE, dtype=np.uint32)
    dq.flat[:50] = pixel.DO_NOT_USE
    return data, dq, None, 5.0, 0.0, 0.5


def _case_non_do_not_use_flags_only():
    data = np.full(SHAPE, 10.0, dtype=np.float32)
    dq = np.full(SHAPE, pixel.SATURATED | pixel.HOT, dtype=np.uint32)
    return data, dq, None, 10.0, 0.0, 1.0


def _case_nan_pixels_excluded_from_statistics():
    data = np.full(SHAPE, 8.0, dtype=np.float32)
    data[0, 0] = np.nan
    dq = np.zeros(SHAPE, dtype=np.uint32)
    expected_fraction = (data.size - 1) / data.size
    return data, dq, None, 8.0, 0.0, expected_fraction


def _case_nonpositive_err_pixels_excluded_from_statistics():
    data = np.full(SHAPE, 7.0, dtype=np.float32)
    dq = np.zeros(SHAPE, dtype=np.uint32)
    err = np.ones(SHAPE, dtype=np.float32)
    err.flat[:40] = 0.0
    expected_fraction = (data.size - 40) / data.size
    return data, dq, err, 7.0, 0.0, expected_fraction


@pytest.mark.parametrize(
    "model_class",
    [ImageModel, MosaicModel],
    ids=["image_model", "mosaic_model"],
)
@pytest.mark.parametrize(
    "case_factory",
    [
        pytest.param(_case_all_pixels_good, id="all_pixels_good"),
        pytest.param(_case_all_pixels_do_not_use, id="all_pixels_do_not_use"),
        pytest.param(_case_half_pixels_do_not_use, id="half_pixels_do_not_use"),
        pytest.param(_case_non_do_not_use_flags_only, id="non_do_not_use_flags_only"),
        pytest.param(
            _case_nan_pixels_excluded_from_statistics,
            id="nan_pixels_excluded_from_statistics",
        ),
        pytest.param(
            _case_nonpositive_err_pixels_excluded_from_statistics,
            id="nonpositive_err_pixels_excluded_from_statistics",
        ),
    ],
)
def test_populate_statistics(model_class, case_factory):
    """Populate statistics with different pixel-validity conditions."""
    data, dq, err, expected_median, expected_rms, expected_frac = case_factory()
    model = _create_model(model_class, data=data, dq=dq, err=err)

    populate_statistics(model)

    _assert_float_or_nan(model.meta.statistics.image_median, expected_median)
    _assert_float_or_nan(model.meta.statistics.image_rms, expected_rms)
    assert model.meta.statistics.good_pixel_fraction == pytest.approx(
        expected_frac, abs=1e-6
    )
    _assert_model_specific_placeholders(model)


@pytest.mark.parametrize(
    "model_class",
    [ImageModel, MosaicModel],
    ids=["image_model", "mosaic_model"],
)
def test_populate_statistics_warns_when_no_good_pixels(model_class, caplog):
    """No-good-pixel path logs a warning and leaves stats at defaults."""
    data = np.full(SHAPE, 12.0, dtype=np.float32)
    dq = np.zeros(SHAPE, dtype=np.uint32)
    err = np.zeros(SHAPE, dtype=np.float32)
    model = _create_model(model_class, data=data, dq=dq, err=err)

    with caplog.at_level("WARNING", logger="romancal.lib.statistics"):
        populate_statistics(model)

    assert "No good pixels found for statistics calculation." in caplog.text
    assert np.isnan(model.meta.statistics.image_median)
    assert np.isnan(model.meta.statistics.image_rms)
    assert model.meta.statistics.good_pixel_fraction == 0.0
    _assert_model_specific_placeholders(model)


@pytest.mark.parametrize(
    "model_class",
    [ImageModel, MosaicModel],
    ids=["image_model", "mosaic_model"],
)
def test_statistics_graceful_exit_no_data(model_class):
    """Ensure meta.statistics is created with defaults if data is None."""
    model = model_class.create_minimal()
    model.data = None

    populate_statistics(model)

    assert hasattr(model.meta, "statistics")
    assert np.isnan(model.meta.statistics.image_median)
    assert np.isnan(model.meta.statistics.image_rms)
    assert model.meta.statistics.good_pixel_fraction == 0.0
    _assert_model_specific_placeholders(model)
