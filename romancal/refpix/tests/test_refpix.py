import numpy as np
from numpy.testing import assert_allclose

from romancal.refpix.data import Aligned, ChannelFFT, Coefficients, Standard
from romancal.refpix.refpix import (
    amp33_cosine_interpolation,
    amp33_fft_interpolation,
    correction,
    remove_linear_trends,
    remove_offset,
)

from .conftest import Dims


class TestOffset:
    def test_function(self, standard_data: Standard):
        new = remove_offset(standard_data)

        # check offset is set
        assert new.offset is not None
        assert new.offset.shape == (Dims.N_ROWS, Dims.N_COLS)
        assert (new.offset != 0).all()

        # check data is updated
        assert (new.data != standard_data.data).all()
        assert (new.amp33 != standard_data.amp33).all()

        # add offset back (error accumulates due to floating point arithmetic)
        reset = Standard.from_combined(new.combined_data + new.offset)
        assert_allclose(reset.data, standard_data.data, atol=1e-5)
        assert_allclose(reset.amp33, standard_data.amp33, atol=1e-5)

    def test_regression(self, standard_data: Standard):
        from . import reference_utils

        # Run regression utility (modifies in place)
        regression = standard_data.combined_data
        _, b = reference_utils.remove_linear_trends(regression, True)

        test = remove_offset(standard_data)

        # Regression test
        assert (test.offset == b).all()
        assert (test.combined_data == regression).all()


class TestRemoveLinearTrends:
    def test_regression(self, aligned_data: Aligned):
        from . import reference_utils

        # Run regression utility
        regression = aligned_data.combined_data
        reference_utils.remove_linear_trends_per_frame(regression, False, True)

        new = remove_linear_trends(aligned_data)
        assert (new.combined_data == regression).all()


class TestAmp33CosineInterpolation:
    def test_regression(self, aligned_data: Aligned):
        from . import reference_utils

        # Run regression utility
        regression = aligned_data.combined_data
        reference_utils.cos_interp_reference(regression, regression.shape[1])
        assert (
            regression[-1, :, :, :] != aligned_data.combined_data[-1, :, :, :]
        ).any()
        assert (
            regression[:-1, :, :, :] == aligned_data.combined_data[:-1, :, :, :]
        ).all()

        new = amp33_cosine_interpolation(aligned_data)

        assert (new.amp33 == regression[-1, :, :, :]).all()
        assert (new.combined_data == regression).all()


class TestAmp33FftInterpolation:
    def test_regression(self, aligned_data: Aligned):
        """
        Run fft interpolation regression test

        NOTE:
            The reference code assumes the data will be changed in-place for all
            its major operations.However, the reference code violates this assumption
            for the Amp33 FFT interpolation step. It does make an in-place change to
            a sub-array, `dataReferenceChan_FramesFlat`, of the main data array, but
            this sub-array does not map to the original data array, `dataUniformTime`.
            The sub-array is not used in the rest of the reference code, which I believe
            is a bug.

            To combat this apparent mistake, I reshape the
            `dataReferenceChan_FramesFlat`
            array in my wrapper of the reference code and then return that array as a
            function output.

            The romancal code does not make this mistake because it does not work
            via in-place changes.
        """
        from . import reference_utils

        regression = aligned_data.combined_data
        amp33_regression = reference_utils.fft_interp_amp33(
            regression, regression.shape[1]
        )

        # Show this sub array does get inplace changes (returned as function output)
        assert (
            amp33_regression[:, :, :] != aligned_data.combined_data[-1, :, :, :]
        ).any()

        # Demonstration of the issue with the reference code, if an in-place change
        # is actually made, then the first assert below would fail.
        assert (regression[-1, :, :, :] == aligned_data.amp33).all()
        assert (regression[:-1, :, :, :] == aligned_data.data).all()

        new = amp33_fft_interpolation(aligned_data, 3)

        assert (amp33_regression == new.amp33).all()


class TestComputeCorrection:
    def test_correction_regression(
        self, channel_fft: ChannelFFT, coefficients: Coefficients
    ):
        from . import reference_utils

        regression = reference_utils.compute_correction(
            coefficients.alpha,
            coefficients.gamma,
            coefficients.zeta,
            Dims.N_FRAMES,
            channel_fft.left,
            channel_fft.right,
            channel_fft.amp33,
        ).real

        correct = np.array(list(correction(channel_fft, coefficients)))
        assert (correct[:, :, :] == regression[:, :, :]).all()

    def test_apply_regression(self, aligned_data: Aligned, coefficients: Coefficients):
        from . import reference_utils

        left = aligned_data.left.copy()
        right = aligned_data.right.copy()
        amp33 = aligned_data.amp33.copy()

        reference_utils.forward_fft(Dims.N_FRAMES, aligned_data.combined_data)

        assert (left == aligned_data.left).all()
        assert (right == aligned_data.right).all()
        assert (amp33 == aligned_data.amp33).all()

        # channel_fft = ChannelFFT.from_aligned(aligned_data)

        # assert (channel_fft.left == regression[0]).all()
        # assert (channel_fft.right == regression[1]).all()
        # assert (channel_fft.amp33 == regression[2]).all()

        # channel_fft2 = ChannelFFT.from_aligned(aligned_data)

        # assert (channel_fft2.left == regression[0]).all()
        # assert (channel_fft2.right == regression[1]).all()
        # assert (channel_fft2.amp33 == regression[2]).all()

        # regression = reference_utils.apply_correction(
        #     Dims.N_FRAMES,
        #     aligned_data.combined_data,
        #     coefficients.alpha,
        #     coefficients.gamma,
        #     coefficients.zeta,
        # )

        # new = apply_correction(aligned_data, coefficients)
        # assert (new.left == regression[0]).all()
        # assert (new.right == regression[1]).all()
        # assert (new.amp33 == regression[2]).all()
        # assert (new == regression).all()

        # assert (new.combined_data[:, :, :, : - Width.PAD] == regression).all()
