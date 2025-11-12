"""
Unit tests for the source catalog utils module.
"""

import numpy as np
import pytest
from roman_datamodels.datamodels import ImageModel, MosaicModel

from romancal.source_catalog.utils import copy_model_arrays

rng = np.random.default_rng(12345)


class TestCopyModelArrays:
    """
    Tests for the copy_model_arrays function.
    """

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_data_copied(self, model_class):
        """
        Test that data array is copied independently.
        """
        model = model_class.create_fake_data(shape=(50, 50))
        model.data[:] = rng.normal(0, 1, size=(50, 50))

        copied = copy_model_arrays(model)

        # Verify data is copied (different memory addresses)
        assert id(model.data) != id(copied.data)
        # Verify data values are equal
        assert np.array_equal(model.data, copied.data)

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_err_copied(self, model_class):
        """
        Test that err array is copied independently.
        """
        model = model_class.create_fake_data(shape=(50, 50))
        model.err = np.ones((50, 50)) * 0.5

        copied = copy_model_arrays(model)

        # Verify err is copied (different memory addresses)
        assert id(model.err) != id(copied.err)
        # Verify err values are equal
        assert np.array_equal(model.err, copied.err)

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_meta_preserved(self, model_class):
        """
        Test that metadata is preserved.
        """
        model = model_class.create_fake_data(shape=(50, 50))

        copied = copy_model_arrays(model)

        # Verify meta has same values (note: datamodels may copy the dict)
        assert model.meta.model_type == copied.meta.model_type
        assert model.meta.telescope == copied.meta.telescope

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_data_modification_independent(self, model_class):
        """
        Test that modifying copied data doesn't affect original.
        """
        model = model_class.create_fake_data(shape=(50, 50))
        model.data[:] = rng.normal(0, 1, size=(50, 50))
        original_value = model.data[0, 0]

        copied = copy_model_arrays(model)
        copied.data[0, 0] = 999.0

        # Original should be unchanged
        assert model.data[0, 0] == original_value
        # Copied should have new value
        assert copied.data[0, 0] == 999.0

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_err_modification_independent(self, model_class):
        """
        Test that modifying copied err doesn't affect original.
        """
        model = model_class.create_fake_data(shape=(50, 50))
        model.err = np.ones((50, 50)) * 0.5
        original_err_value = model.err[10, 10]

        copied = copy_model_arrays(model)
        copied.err[10, 10] = 123.456

        # Original should be unchanged
        assert model.err[10, 10] == original_err_value
        # Copied should have new value
        assert copied.err[10, 10] == 123.456

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_returns_correct_type(self, model_class):
        """
        Test that function returns the same model type.
        """
        model = model_class.create_fake_data(shape=(50, 50))
        copied = copy_model_arrays(model)
        assert isinstance(copied, model_class)

    @pytest.mark.parametrize("model_class", [ImageModel, MosaicModel])
    def test_preserves_data_shape(self, model_class):
        """
        Test that copied model preserves data shape.
        """
        shape = (100, 75)
        model = model_class.create_fake_data(shape=shape)
        copied = copy_model_arrays(model)
        assert copied.data.shape == shape

    def test_image_model_dq_shared(self):
        """
        Test that ImageModel dq array is shared (not copied).
        """
        model = ImageModel.create_fake_data(shape=(50, 50))

        copied = copy_model_arrays(model)

        # Verify dq is shared (same memory address)
        assert id(model.dq) == id(copied.dq)

    def test_mosaic_model_weight_shared(self):
        """
        Test that MosaicModel weight array is shared (not copied).
        """
        model = MosaicModel.create_fake_data(shape=(50, 50))
        model.weight = np.ones((50, 50))

        copied = copy_model_arrays(model)

        # Verify weight is shared (same memory address)
        assert id(model.weight) == id(copied.weight)

    def test_invalid_model_type(self):
        """
        Test that invalid model type raises TypeError.
        """

        class FakeModel:
            pass

        expected_msg = "model must be an ImageModel or MosaicModel"
        with pytest.raises(TypeError, match=expected_msg):
            copy_model_arrays(FakeModel())
