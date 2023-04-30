import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose
from roman_datamodels.maker_utils import mk_ramp

from romancal.refpix.refpix import Arrangement, RefPixData, Width

RNG = np.random.default_rng(42)
N_FRAMES = 8
N_ROWS = 4096
N_COLS = N_ROWS + Width.CHANNEL
N_CHAN = N_COLS // Width.CHANNEL
N_DETECT_CHAN = N_CHAN - 1


@pytest.fixture(scope="module")
def ramp_model():

    datamodel = mk_ramp()
    assert datamodel.data.shape == (N_FRAMES, N_ROWS, N_ROWS)
    assert datamodel.amp33.shape == (N_FRAMES, N_ROWS, Width.CHANNEL)

    data = RNG.uniform(1, 100, size=(N_FRAMES, N_ROWS, N_ROWS)).astype(
        datamodel.data.dtype
    )
    amp33 = RNG.uniform(1, 100, size=(N_FRAMES, N_ROWS, Width.CHANNEL)).astype(
        datamodel.amp33.dtype
    )

    datamodel.data = u.Quantity(
        data, unit=datamodel.data.unit, dtype=datamodel.data.dtype
    )
    datamodel.border_ref_pix_left = datamodel.data[:, :, : Width.REF]
    datamodel.border_ref_pix_right = datamodel.data[:, :, -Width.REF :]

    datamodel.amp33 = u.Quantity(
        amp33, unit=datamodel.amp33.unit, dtype=datamodel.amp33.dtype
    )

    return datamodel


def test_construct_data(ramp_model):
    datamodel = RefPixData.from_datamodel(ramp_model)

    assert (datamodel.data == ramp_model.data.value).all()
    assert (datamodel.amp33 == ramp_model.amp33.value).all()
    assert (datamodel.left == ramp_model.border_ref_pix_left.value).all()
    assert (datamodel.right == ramp_model.border_ref_pix_right.value).all()
    assert datamodel.offset is None

    assert datamodel.data.shape == (N_FRAMES, N_ROWS, N_ROWS)
    assert datamodel.amp33.shape == (N_FRAMES, N_ROWS, Width.CHANNEL)
    assert datamodel.left.shape == (N_FRAMES, N_ROWS, Width.REF)
    assert datamodel.right.shape == (N_FRAMES, N_ROWS, Width.REF)

    assert (datamodel.left == datamodel.data[:, :, : Width.REF]).all()
    assert (datamodel.right == datamodel.data[:, :, -Width.REF :]).all()

    assert datamodel.arrangement == Arrangement.STANDARD


@pytest.fixture()
def ref_data(ramp_model):
    return RefPixData.from_datamodel(ramp_model)


def test_combine_data(ref_data):
    data = ref_data.combine_data

    assert data.shape == (N_FRAMES, N_ROWS, N_COLS)
    assert (data[:, :, :N_ROWS] == ref_data.data).all()
    assert (data[:, :, N_ROWS:] == ref_data.amp33).all()

    # Check function in split mode
    new = RefPixData.from_split(ref_data.split_channels)
    assert (new.combine_data == data).all()


def test_from_combined_data(ref_data):
    pert = 37

    data = ref_data.data + pert
    amp33 = ref_data.amp33 + pert

    # apply the perturbation to combined data then unpack
    new = RefPixData.from_combined_data(ref_data.combine_data + pert)
    assert (new.data == data).all()
    assert (new.amp33 == amp33).all()
    assert new.arrangement == Arrangement.STANDARD


def test_remove_offset(ref_data):
    new = ref_data.remove_offset()

    # check offset is set
    assert new.offset is not None
    assert new.offset.shape == (N_ROWS, N_COLS)
    assert (new.offset != 0).all()

    # check data is updated
    assert (new.data != ref_data.data).all()
    assert (new.amp33 != ref_data.amp33).all()

    # add offset back (error accumulates due to floating point arithmetic)
    reset = RefPixData.from_combined_data(new.combine_data + new.offset)
    assert_allclose(reset.data, ref_data.data, atol=1e-5)
    assert_allclose(reset.amp33, ref_data.amp33, atol=1e-5)


def test_regress_remove_offset(ref_data):
    from . import reference_utils

    # Run regression utility (modifies in place)
    data_regress = ref_data.combine_data
    _, b = reference_utils.remove_linear_trends(data_regress, True)

    test = ref_data.remove_offset()

    # Regression test
    assert (test.offset == b).all()
    assert (test.combine_data == data_regress).all()


def test_split_channels(ref_data):
    split = ref_data.split_channels

    assert split.shape == (N_CHAN, N_FRAMES, N_ROWS, Width.CHANNEL)
    # check amp33 is preserved
    assert (split[-1, :, :, :] == ref_data.amp33).all()

    # check data is split correctly
    for idx in range(0, N_DETECT_CHAN):
        assert (
            split[idx, :, :, 0 : Width.CHANNEL]
            == ref_data.data[:, :, idx * Width.CHANNEL : (idx + 1) * Width.CHANNEL]
        ).all()

    # Check function in split mode
    new = RefPixData.from_split(split)
    assert (new.split_channels == split).all()


def test_aligned_channels(ref_data):
    data = ref_data.split_channels

    aligned = ref_data.aligned_channels
    # Check shape includes padding
    assert aligned.shape == (N_CHAN, N_FRAMES, N_ROWS, Width.CHANNEL + Width.PAD)

    # check that it is Padded Correctly
    assert (aligned[3:, :, :, -Width.PAD :] == 0).all()

    # Check data is reversed and split correctly
    # Even channels are the same
    assert (aligned[::2, :, :, : Width.CHANNEL] == data[::2, :, :, :]).all()
    # Odd channels are flipped
    assert (aligned[1::2, :, :, : Width.CHANNEL] == data[1::2, :, :, ::-1]).all()


def test_regress_aligned_channels(ref_data):
    from . import reference_utils

    # Run regression utility
    split_regress = reference_utils.aligned_channels(ref_data.combine_data)

    # Regression test
    split = ref_data.aligned_channels
    assert (split_regress == split).all()


def test_from_split(ref_data):
    offset = np.array([1, 2, 3])
    data = ref_data.aligned_channels

    new = RefPixData.from_split(data, offset)
    assert (new.data == data[:-1, :, :, :]).all()
    assert (new.amp33 == data[-1, :, :, :]).all()
    assert (new.offset == offset).all()
    assert new.arrangement == Arrangement.SPLIT


def test_regress_remove_linear_trends(ref_data):
    from . import reference_utils

    # Run regression utility
    regress = ref_data.aligned_channels
    reference_utils.remove_linear_trends_per_frame(regress, False, True)

    new = ref_data.remove_linear_trends()
    assert_allclose(new, regress)
