import numpy as np
import pytest
from astropy import units as u
from roman_datamodels.maker_utils import mk_ramp

from romancal.refpix.data import Aligned, Standard, Width

RNG = np.random.default_rng(42)
N_FRAMES = 8
N_ROWS = 4096
N_COLS = N_ROWS + Width.CHANNEL
N_CHAN = N_COLS // Width.CHANNEL
N_DETECT_CHAN = N_CHAN - 1
N_ALIGN_COLS = Width.CHANNEL + Width.PAD


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
    datamodel = Standard.from_datamodel(ramp_model)

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


@pytest.fixture(scope="module")
def standard_data(ramp_model) -> Standard:
    return Standard.from_datamodel(ramp_model)


@pytest.fixture(scope="module")
def aligned_data(standard_data) -> Aligned:
    return Aligned.from_standard(standard_data)


class TestStandard:
    def test_combined_data(self, standard_data):
        data = standard_data.combined_data

        assert data.shape == (N_FRAMES, N_ROWS, N_COLS)
        assert (data[:, :, :N_ROWS] == standard_data.data).all()
        assert (data[:, :, N_ROWS:] == standard_data.amp33).all()

    def test_from_combined(self, standard_data):
        data = standard_data.combined_data

        new = Standard.from_combined(data)
        assert (new.data == standard_data.data).all()
        assert (new.amp33 == standard_data.amp33).all()

    def test_aligned_data(self, standard_data):
        data = standard_data.aligned_data

        # Check the shape is correct (including padding)
        assert data.shape == (N_CHAN, N_FRAMES, N_ROWS, N_ALIGN_COLS)

        # Check amp33 is preserved
        assert (data[-1, :, :, : Width.CHANNEL] == standard_data.amp33).all()

        # Check that it is Padded Correctly
        assert (data[:, :, :, -Width.PAD :] == 0).all()

        for idx in range(0, N_CHAN):
            truth = standard_data.combined_data[
                :, :, idx * Width.CHANNEL : (idx + 1) * Width.CHANNEL
            ]
            if idx % 2 == 0:
                # Even channels are not flipped
                assert (data[idx, :, :, : Width.CHANNEL] == truth).all()
            else:
                # Odd channels are flipped on last axis
                assert (data[idx, :, :, : Width.CHANNEL] == truth[:, :, ::-1]).all()

    def test_from_aligned(self, aligned_data):
        new = Standard.from_aligned(aligned_data)

        # Check the shapes are correct
        assert new.data.shape == (N_FRAMES, N_ROWS, N_ROWS)
        assert new.amp33.shape == (N_FRAMES, N_ROWS, Width.CHANNEL)

        # Check that the data is transformed back correctly
        for idx in range(0, N_DETECT_CHAN):
            test_data = new.data[:, :, idx * Width.CHANNEL : (idx + 1) * Width.CHANNEL]
            if idx % 2 == 0:
                # Even channels are not flipped
                assert (
                    test_data == aligned_data.data[idx, :, :, : Width.CHANNEL]
                ).all()
            else:
                # Odd channels are flipped on last axis
                assert (
                    test_data[:, :, ::-1]
                    == aligned_data.data[idx, :, :, : Width.CHANNEL]
                ).all()

        # Check that the amp33 is transformed correctly
        assert (new.amp33 == aligned_data.amp33[:, :, : Width.CHANNEL]).all()


class TestAligned:
    def test_combined_data(self, aligned_data):
        data = aligned_data.combined_data

        assert data.shape == (N_CHAN, N_FRAMES, N_ROWS, N_ALIGN_COLS)
        assert (data[:-1, :, :, :] == aligned_data.data).all()
        assert (data[-1, :, :, :] == aligned_data.amp33).all()

    def test_from_combined(self, aligned_data):
        data = aligned_data.combined_data

        new = Aligned.from_combined(data)
        assert (new.data == aligned_data.data).all()
        assert (new.amp33 == aligned_data.amp33).all()

    def test_standard_data(self, aligned_data):
        data = aligned_data.standard_data

        # Check shape is correct (and padding is removed)
        assert data.shape == (N_FRAMES, N_ROWS, N_COLS)

        for idx in range(0, N_CHAN):
            test_data = data[:, :, idx * Width.CHANNEL : (idx + 1) * Width.CHANNEL]
            if idx % 2 == 0:
                # Even channels are not flipped
                assert (
                    test_data == aligned_data.combined_data[idx, :, :, : Width.CHANNEL]
                ).all()
            else:
                # Odd channels are flipped on last axis
                assert (
                    test_data[:, :, ::-1]
                    == aligned_data.combined_data[idx, :, :, : Width.CHANNEL]
                ).all()

    def test_from_standard(self, standard_data):
        new = Aligned.from_standard(standard_data)

        # Check the shapes are correct
        assert new.data.shape == (N_DETECT_CHAN, N_FRAMES, N_ROWS, N_ALIGN_COLS)
        assert new.amp33.shape == (N_FRAMES, N_ROWS, N_ALIGN_COLS)

        # Check that it is Padded Correctly
        assert (new.data[:, :, :, -Width.PAD :] == 0).all()
        assert (new.amp33[:, :, -Width.PAD :] == 0).all()

        # Check data is transformed correctly
        for idx in range(0, N_DETECT_CHAN):
            truth = standard_data.data[
                :, :, idx * Width.CHANNEL : (idx + 1) * Width.CHANNEL
            ]
            if idx % 2 == 0:
                # Even channels are not flipped
                assert (new.data[idx, :, :, : Width.CHANNEL] == truth).all()
            else:
                # Odd channels are flipped on last axis
                assert (new.data[idx, :, :, : Width.CHANNEL] == truth[:, :, ::-1]).all()

        # Check amp33 is transformed correctly
        assert (new.amp33[:, :, : Width.CHANNEL] == standard_data.amp33).all()
