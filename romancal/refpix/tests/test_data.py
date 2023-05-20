import numpy as np

from romancal.refpix.data import Aligned, ChannelFFT, Standard, Width, channel_fft

from .conftest import Dims


def test_construct_data(ramp_model):
    datamodel = Standard.from_datamodel(ramp_model)

    assert (datamodel.data == ramp_model.data.value).all()
    assert (datamodel.amp33 == ramp_model.amp33.value).all()
    assert (datamodel.left == ramp_model.border_ref_pix_left.value).all()
    assert (datamodel.right == ramp_model.border_ref_pix_right.value).all()
    assert datamodel.offset is None

    assert datamodel.data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ROWS)
    assert datamodel.amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.CHANNEL)
    assert datamodel.left.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.REF)
    assert datamodel.right.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.REF)

    assert (datamodel.left == datamodel.data[:, :, : Width.REF]).all()
    assert (datamodel.right == datamodel.data[:, :, -Width.REF :]).all()


class TestStandard:
    def test_combined_data(self, standard_data: Standard):
        data = standard_data.combined_data

        assert data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)
        assert (data[:, :, : Dims.N_ROWS] == standard_data.data).all()
        assert (data[:, :, Dims.N_ROWS :] == standard_data.amp33).all()

    def test_from_combined(self, standard_data: Standard):
        data = standard_data.combined_data

        new = Standard.from_combined(data)
        assert (new.data == standard_data.data).all()
        assert (new.amp33 == standard_data.amp33).all()

    def test_aligned_data(self, standard_data):
        data = standard_data.aligned_data

        # Check the shape is correct (including padding)
        assert data.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.N_ALIGN_COLS,
        )

        # Check amp33 is preserved
        assert (data[-1, :, :, : Width.CHANNEL] == standard_data.amp33).all()

        # Check that it is Padded Correctly
        assert (data[:, :, :, -Width.PAD :] == 0).all()

        for idx in range(0, Dims.N_CHAN):
            truth = standard_data.combined_data[
                :, :, idx * Width.CHANNEL : (idx + 1) * Width.CHANNEL
            ]
            if idx % 2 == 0:
                # Even channels are not flipped
                assert (data[idx, :, :, : Width.CHANNEL] == truth).all()
            else:
                # Odd channels are flipped on last axis
                assert (data[idx, :, :, : Width.CHANNEL] == truth[:, :, ::-1]).all()

    def test_from_aligned(self, aligned_data: Aligned):
        new = Standard.from_aligned(aligned_data)

        # Check the shapes are correct
        assert new.data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ROWS)
        assert new.amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.CHANNEL)

        # Check that the data is transformed back correctly
        for idx in range(0, Dims.N_DETECT_CHAN):
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

    def test_roundtrip(self, standard_data: Standard):
        # Insurance that the data in standard_data is not modified
        data = standard_data.data.copy()
        amp33 = standard_data.amp33.copy()

        # Convert to aligned and then back
        aligned = Aligned.from_standard(standard_data)
        new = Standard.from_aligned(aligned)

        assert (new.data == standard_data.data).all()
        assert (new.amp33 == standard_data.amp33).all()
        assert (new.data == data).all()
        assert (new.amp33 == amp33).all()

    def test_regression(self, standard_data: Standard):
        """Regression test to the reference data alignment"""
        from . import reference_utils

        # Run regression utility
        regression = reference_utils.aligned_channels(standard_data.combined_data)

        # Regression test
        assert (regression == standard_data.aligned_data).all()

    def test_left(self, standard_data: Standard):
        left = standard_data.left

        assert left.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.REF)
        assert (left == standard_data.data[:, :, : Width.REF]).all()

    def test_right(self, standard_data: Standard):
        right = standard_data.right

        assert right.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.REF)
        assert (right == standard_data.data[:, :, -Width.REF :]).all()


class TestAligned:
    def test_combined_data(self, aligned_data: Aligned):
        data = aligned_data.combined_data

        assert data.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.N_ALIGN_COLS,
        )
        assert (data[:-1, :, :, :] == aligned_data.data).all()
        assert (data[-1, :, :, :] == aligned_data.amp33).all()

    def test_from_combined(self, aligned_data: Aligned):
        data = aligned_data.combined_data

        new = Aligned.from_combined(data)
        assert (new.data == aligned_data.data).all()
        assert (new.amp33 == aligned_data.amp33).all()

    def test_standard_data(self, aligned_data: Aligned):
        data = aligned_data.standard_data

        # Check shape is correct (and padding is removed)
        assert data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)

        for idx in range(0, Dims.N_CHAN):
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

    def test_from_standard(self, standard_data: Standard):
        new = Aligned.from_standard(standard_data)

        # Check the shapes are correct
        assert new.data.shape == (
            Dims.N_DETECT_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.N_ALIGN_COLS,
        )
        assert new.amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ALIGN_COLS)

        # Check that it is Padded Correctly
        assert (new.data[:, :, :, -Width.PAD :] == 0).all()
        assert (new.amp33[:, :, -Width.PAD :] == 0).all()

        # Check data is transformed correctly
        for idx in range(0, Dims.N_DETECT_CHAN):
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

    def test_roundtrip(self, aligned_data: Aligned):
        # Insurance that the data in standard_data is not modified
        data = aligned_data.data.copy()
        amp33 = aligned_data.amp33.copy()

        # Convert to aligned and then back
        standard = Standard.from_aligned(aligned_data)
        new = Aligned.from_standard(standard)

        assert (new.data == aligned_data.data).all()
        assert (new.amp33 == aligned_data.amp33).all()
        assert (new.data == data).all()
        assert (new.amp33 == amp33).all()

    def test_regression(self, aligned_data: Aligned):
        """Regression test to the reference data alignment"""
        from . import reference_utils

        # Run regression utility
        regression = reference_utils.aligned_channels(aligned_data.standard_data)

        # Regression test
        assert (regression == aligned_data.combined_data).all()

    def test_left(self, aligned_data: Aligned):
        left = aligned_data.left

        # Check values
        assert left.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ALIGN_COLS)
        assert (
            left[:, :, : Width.REF] == aligned_data.data[0, :, :, : Width.REF]
        ).all()
        assert (left[:, :, Width.REF :] == 0).all()

        # Check against the standard form of the data
        standard = Standard.from_aligned(aligned_data)
        assert (left[:, :, : Width.REF] == standard.left).all()

    def test_left_regression(self, aligned_data: Aligned):
        from . import reference_utils

        # Insurance against in-place modifications
        left = aligned_data.left.copy()

        regression = reference_utils.left(aligned_data.combined_data)
        assert (regression == aligned_data.left).all()
        assert (regression == left).all()

    def test_right(self, aligned_data: Aligned):
        right = aligned_data.right

        # Check values
        assert right.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_ALIGN_COLS)
        assert (
            right[:, :, -Width.REF :] == aligned_data.data[-1, :, :, -Width.REF :]
        ).all()
        assert (right[:, :, : -Width.REF] == 0).all()

        # Check against the standard form of the data
        standard = Standard.from_aligned(aligned_data)
        assert (right[:, :, -Width.REF :] == standard.right).all()

    def test_right_regression(self, aligned_data: Aligned):
        from . import reference_utils

        # # Insurance against in-place modifications
        # right = aligned_data.right.copy()
        # Demonstrate the correct implementation
        regression = reference_utils.right(aligned_data.combined_data)
        assert (regression == aligned_data.right).all()
        # assert (regression == right).all()

    def test_bug_right(self, standard_data: Aligned):
        from . import reference_utils

        data = standard_data.combined_data.copy()

        # Original Right:
        right = np.copy(data[:, :, -Width.REF - Width.CHANNEL : -Width.CHANNEL])
        assert right.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.REF)
        assert (right == standard_data.right).all()

        # align data:
        aligned = reference_utils.aligned_channels(data)
        assert aligned.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.N_ALIGN_COLS,
        )
        right_channel_padded = np.copy(aligned[-2, :, :, :])
        assert right_channel_padded.shape == (
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.N_ALIGN_COLS,
        )
        assert (right_channel_padded[:, :, -Width.PAD :] == 0).all()
        right_channel = np.copy(right_channel_padded[:, :, : -Width.PAD])
        assert right_channel.shape == (Dims.N_FRAMES, Dims.N_ROWS, Width.CHANNEL)
        right_channel_flip = np.copy(right_channel[:, :, ::-1])
        assert (right_channel_flip[:, :, -4:] == right).all()

        assert (right_channel[:, :, :4].flip() == right).all()

        # assert (right_channel[:, :, -Width.REF - Width.PAD:-Width.PAD] == right).all()

        # aligned_right = reference_utils.right(aligned)
        # assert (aligned_right[:, :, :4] == right).all()


class TestChannelFFT:
    def test_left_regression(self, aligned_data: Aligned):
        from . import reference_utils

        left = channel_fft(aligned_data.left, normalize=True)
        regression = reference_utils.left_fft(Dims.N_FRAMES, aligned_data.combined_data)
        assert (left == regression).all()

    def test_right_regression(self, aligned_data: Aligned):
        from . import reference_utils

        right = channel_fft(aligned_data.right, normalize=True)
        regression = reference_utils.right_fft(
            Dims.N_FRAMES, aligned_data.combined_data
        )
        assert (right == regression).all()

    def test_amp33_regression(self, aligned_data: Aligned):
        from . import reference_utils

        amp33 = channel_fft(aligned_data.amp33)
        regression = reference_utils.amp33_fft(
            Dims.N_FRAMES, aligned_data.combined_data
        )
        assert (amp33 == regression).all()

    def test_from_aligned(self, aligned_data: Aligned):
        from . import reference_utils

        # left = aligned_data.left.copy()
        # right = aligned_data.right.copy()
        # amp33 = aligned_data.amp33.copy()
        aligned_data.data.copy()

        regression = reference_utils.forward_fft(
            Dims.N_FRAMES, aligned_data.combined_data
        )
        # assert (left == aligned_data.left).all()
        # assert (right == aligned_data.right).all()
        # assert (amp33 == aligned_data.amp33).all()
        regression[0].copy()
        regression[1].copy()
        regression[2].copy()
        channel_fft = ChannelFFT.from_aligned(aligned_data)
        # assert (r_left == regression[0]).all()
        # assert (r_right == regression[1]).all()
        # assert (r_amp33 == regression[2]).all()
        # assert (left == aligned_data.left).all()
        # assert (right == aligned_data.right).all()
        # assert (amp33 == aligned_data.amp33).all()

        assert (channel_fft.left == regression[0]).all()
        assert (channel_fft.right == regression[1]).all()
        assert (channel_fft.amp33 == regression[2]).all()
