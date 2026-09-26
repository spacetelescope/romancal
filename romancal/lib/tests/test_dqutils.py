"""Tests for the data quality propagation helpers."""

import numpy as np
import pytest
import roman_datamodels.datamodels as rdm

from romancal.lib.dqutils import DQ2_DTYPE, dq_names, ensure_dq2, update_dq

SHAPE = (16, 16)
UNTRIMMED = (24, 24)


def make_ramp():
    model = rdm.RampModel.create_fake_data(shape=(3, *SHAPE))
    model.pixeldq = np.zeros(SHAPE, dtype=np.uint32)
    model.pop("pixeldq2", None)
    return model


def make_image():
    model = rdm.ImageModel.create_fake_data(shape=SHAPE)
    model.dq = np.zeros(SHAPE, dtype=np.uint32)
    model.pop("dq2", None)
    return model


def make_mask(shape=SHAPE, dq=1, dq2=None):
    model = rdm.MaskRefModel.create_fake_data(shape=shape)
    model.dq = np.full(shape, dq, dtype=np.uint32)
    if dq2 is None:
        model.pop("dq2", None)
    else:
        model["dq2"] = np.full(shape, dq2, dtype=DQ2_DTYPE)
    return model


def test_dq_names_dispatches_on_model_type():
    assert dq_names(make_ramp()) == ("pixeldq", "pixeldq2")
    assert dq_names(make_image()) == ("dq", "dq2")


def test_ensure_dq2_creates_matching_array():
    model = make_image()
    assert "dq2" not in model

    dq2 = ensure_dq2(model)

    assert dq2.shape == model.dq.shape
    assert dq2.dtype == DQ2_DTYPE
    assert not dq2.any()


def test_ensure_dq2_preserves_existing_array():
    model = make_image()
    model["dq2"] = np.full(SHAPE, 7, dtype=DQ2_DTYPE)

    assert ensure_dq2(model) is model.dq2
    assert (model.dq2 == 7).all()


def test_ensure_dq2_ignores_models_without_dq():
    model = rdm.ScienceRawModel.create_fake_data(shape=(3, *SHAPE))

    assert ensure_dq2(model) is None
    assert "dq2" not in model


def test_update_dq_uses_ramp_array_names():
    model = make_ramp()

    update_dq(model, make_mask(dq=5, dq2=9))

    assert (model.pixeldq == 5).all()
    assert (model.pixeldq2 == 9).all()
    assert "dq2" not in model


def test_update_dq_creates_dq2_when_reference_lacks_it():
    model = make_image()

    update_dq(model, make_mask(dq=3))

    assert (model.dq == 3).all()
    assert not model.dq2.any()
    assert model.dq2.dtype == DQ2_DTYPE


def test_update_dq_ors_rather_than_overwrites():
    model = make_image()
    model.dq[:] = 2
    model["dq2"] = np.full(SHAPE, 4, dtype=DQ2_DTYPE)

    update_dq(model, make_mask(dq=1, dq2=8))

    assert (model.dq == 3).all()
    assert (model.dq2 == 12).all()


def test_update_dq_is_idempotent():
    model = make_image()
    reference = make_mask(dq=6, dq2=10)

    update_dq(model, reference)
    expected_dq = model.dq.copy()
    expected_dq2 = model.dq2.copy()
    update_dq(model, reference)

    np.testing.assert_array_equal(model.dq, expected_dq)
    np.testing.assert_array_equal(model.dq2, expected_dq2)


def test_update_dq_preserves_dtypes():
    model = make_image()

    update_dq(model, make_mask(dq=1, dq2=1))

    assert model.dq.dtype == np.uint32
    assert model.dq2.dtype == DQ2_DTYPE


def test_update_dq_ignores_references_without_dq_arrays():
    model = make_image()
    model.dq[:] = 5
    reference = rdm.GainRefModel.create_fake_data(shape=SHAPE)
    assert "dq" not in reference

    update_dq(model, reference)

    assert (model.dq == 5).all()
    assert not model.dq2.any()


def test_update_dq_trims_untrimmed_references():
    """Dark references are full frame while post-ramp-fit science is not."""
    model = make_image()
    reference = make_mask(shape=UNTRIMMED, dq=2, dq2=16)
    # flags that live only in the reference pixel border
    reference.dq[0, 0] = 1 << 10
    reference.dq2[0, 0] = 1 << 11

    update_dq(model, reference)

    assert (model.dq == 2).all()
    assert (model.dq2 == 16).all()


def test_update_dq_does_not_trim_matching_references():
    """Flat references are already trimmed."""
    model = make_image()

    update_dq(model, make_mask(shape=SHAPE, dq=64))

    assert (model.dq == 64).all()


@pytest.mark.parametrize("shape", [(20, 20), (32, 32), (16, 24)])
def test_update_dq_rejects_unexpected_shapes(shape):
    model = make_image()

    with pytest.raises(ValueError, match="matches neither"):
        update_dq(model, make_mask(shape=shape))
