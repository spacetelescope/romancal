"""Tests for the data quality propagation helpers."""

import numpy as np
import pytest
import roman_datamodels.datamodels as rdm

from romancal.lib.dqutils import update_dq

SHAPE = (16, 16)
UNTRIMMED = (24, 24)


def make_image(dq=0, dq2=None):
    model = rdm.ImageModel.create_fake_data(shape=SHAPE)
    model.dq = np.full(SHAPE, dq, dtype=np.uint32)
    if dq2 is None:
        model.pop("dq2", None)
    else:
        model["dq2"] = np.full(SHAPE, dq2, dtype=np.uint32)
    return model


def make_ref(shape=SHAPE, dq=1, dq2=None):
    model = rdm.MaskRefModel.create_fake_data(shape=shape)
    model.dq = np.full(shape, dq, dtype=np.uint32)
    if dq2 is None:
        model.pop("dq2", None)
    else:
        model["dq2"] = np.full(shape, dq2, dtype=np.uint32)
    return model


def test_update_dq_ors_both_arrays():
    """Existing bits are preserved and reference bits are added."""
    model = make_image(dq=2, dq2=4)

    update_dq(model, make_ref(dq=1, dq2=8))

    assert (model.dq == 3).all()
    assert (model.dq2 == 12).all()
    assert model.dq.dtype == np.uint32
    assert model.dq2.dtype == np.uint32


def test_update_dq_uses_ramp_array_names():
    model = rdm.RampModel.create_fake_data(shape=(3, *SHAPE))
    model.pixeldq = np.zeros(SHAPE, dtype=np.uint32)
    model.pop("pixeldq2", None)

    update_dq(model, make_ref(dq=5, dq2=9))

    assert (model.pixeldq == 5).all()
    assert (model.pixeldq2 == 9).all()
    assert "dq2" not in model


@pytest.mark.parametrize("reference", ["mask", "gain"])
def test_update_dq_creates_empty_dq2(reference):
    """dq2 is materialized whether or not the reference supplies one."""
    model = make_image(dq=5)
    if reference == "mask":
        # a reference with dq but no dq2
        reference = make_ref(dq=5)
    else:
        # a reference with neither
        reference = rdm.GainRefModel.create_fake_data(shape=SHAPE)
        assert "dq" not in reference

    update_dq(model, reference)

    assert (model.dq == 5).all()
    assert model.dq2.shape == SHAPE
    assert model.dq2.dtype == np.uint32
    assert not model.dq2.any()


def test_update_dq_ignores_models_without_dq():
    """Models with no pixel-level dq, such as L1, are left alone."""
    model = rdm.ScienceRawModel.create_fake_data(shape=(3, *SHAPE))

    update_dq(model, make_ref())

    assert "dq2" not in model


def test_update_dq_trims_only_untrimmed_references():
    """Dark references are full frame; flat references already are not."""
    model = make_image()
    untrimmed = make_ref(shape=UNTRIMMED, dq=2, dq2=16)
    # bits living only in the reference pixel border, which is dropped
    untrimmed.dq[0, 0] = 1 << 10
    untrimmed.dq2[0, 0] = 1 << 11

    update_dq(model, untrimmed)
    update_dq(model, make_ref(shape=SHAPE, dq=64))

    assert (model.dq == 66).all()
    assert (model.dq2 == 16).all()


@pytest.mark.parametrize("shape", [(20, 20), (32, 32), (16, 24)])
def test_update_dq_rejects_unexpected_shapes(shape):
    model = make_image()

    with pytest.raises(ValueError, match="matches neither"):
        update_dq(model, make_ref(shape=shape))
