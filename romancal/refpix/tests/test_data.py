import numpy as np
from numpy.testing import assert_allclose

from romancal.refpix.data import (
    ChannelView,
    Coefficients,
    Const,
    ReferenceFFT,
    StandardView,
)

from . import reference_utils
from .conftest import RNG, Dims


def test_constants_sanity():
    """
    Sanity check on the constants used and their relation to the test data.
        SMOKE TEST for adjustment to these assumptions.
    """
    assert Const.REF == 4
    assert Const.PAD == 12
    assert Const.CHAN_WIDTH == 128
    assert Const.N_DETECT_CHAN == 32

    # Check the relationship between the constants
    assert Const.N_COLUMNS == 4096
    assert Dims.N_ROWS == 12
    assert Dims.N_COLS == 4224

    # Check that there are the correct total number of channels
    assert Dims.N_CHAN == 33 == Const.N_DETECT_CHAN + 1

    # Check the padded channel width
    assert Dims.PADDED_WIDTH == 140 == Const.CHAN_WIDTH + Const.PAD

    # Check the number of columns in the FFT
    assert Dims.N_FFT_COLS == 841


def test_data(data):
    """
    Sanity that the generated data is the correct shape and type.
    """
    assert data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)
    assert data.dtype == np.float32


class TestStandardView:
    def test_construct_from_datamodel(self, datamodel):
        """
        Test construct from the datamodel passed in
        """
        standard = StandardView.from_datamodel(datamodel)

        # StandardView's data should not be a view to the datamodel
        assert standard.data.base is None

        # Check the relationship between the standard view and the datamodel
        assert (standard.detector == datamodel.data.value).all()
        assert (standard.left == datamodel.border_ref_pix_left.value).all()
        assert (standard.right == datamodel.border_ref_pix_right.value).all()

        # The amp33's dtype changes because it needs to be promoted to match that
        # of the rest of the data
        assert (standard.amp33 == datamodel.amp33.value.astype(np.float32)).all()

    def test_update(self, datamodel):
        standard = StandardView.from_datamodel(datamodel)
        old_data = standard.data.copy()

        # Update the standard view's data
        standard.data = RNG.uniform(
            1, 100, size=(Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)
        ).astype(np.float32)
        assert (standard.data != old_data).any()

        old_detector = datamodel.data.copy()
        old_amp33 = datamodel.amp33.copy()
        old_left = datamodel.border_ref_pix_left.copy()
        old_right = datamodel.border_ref_pix_right.copy()
        # Update the datamodel's data
        new = standard.update(datamodel)

        # Check that new is the datamodel
        assert new is datamodel

        # Check that the datamodel's data has been updated
        assert (new.data != old_detector).any()
        assert (new.amp33 != old_amp33).any()
        assert (new.border_ref_pix_left != old_left).any()
        assert (new.border_ref_pix_right != old_right).any()

        # Check the dtype has been preserved
        assert new.data.dtype == old_detector.dtype
        assert new.amp33.dtype == old_amp33.dtype
        assert new.border_ref_pix_left.dtype == old_left.dtype
        assert new.border_ref_pix_right.dtype == old_right.dtype

        # Check the unit has been preserved
        assert new.data.unit == old_detector.unit
        assert new.amp33.unit == old_amp33.unit
        assert new.border_ref_pix_left.unit == old_left.unit
        assert new.border_ref_pix_right.unit == old_right.unit

        # Check the data has been updated correctly
        assert (new.data.value == standard.detector).all()
        assert (new.border_ref_pix_left.value == standard.left).all()
        assert (new.border_ref_pix_right.value == standard.right).all()

        # The amp33's dtype changes because it needs to be shifted to match the
        # original data's dtype
        assert (new.amp33.value == standard.amp33.astype(old_amp33.dtype)).all()

    def test_create_standard_view(self, data):
        """
        Sanity that the StandardView constructs correctly
        """
        standard = StandardView(data)
        assert standard.data is data
        assert standard.offset is None

    def test_detector(self, data):
        """"""
        standard = StandardView(data)
        detector = standard.detector

        # Check that the detector is the correct shape and type
        assert detector.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.N_COLUMNS)
        assert detector.dtype == np.float32

        # Check that the detector is the correct part of the data
        assert (detector == data[:, :, : Const.N_COLUMNS]).all()

        # Check that the detector is a view of the data
        assert detector.base is data
        assert detector.base is standard.data

    def test_left(self, data):
        standard = StandardView(data)
        left = standard.left

        # Check that the left is the correct shape and type
        assert left.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.REF)
        assert left.dtype == np.float32

        # Check that the left is the correct part of the data
        assert (left == data[:, :, : Const.REF]).all()

        # Check that the left is a view of the data
        assert left.base is data
        assert left.base is standard.data

    def test_right(self, data):
        standard = StandardView(data)
        right = standard.right

        # Check that the right is the correct shape and type
        assert right.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.REF)
        assert right.dtype == np.float32

        # Check that the right is the correct part of the data
        assert (
            right == data[:, :, -Const.REF - Const.CHAN_WIDTH : -Const.CHAN_WIDTH]
        ).all()

        # Check that the right is a view of the data
        assert right.base is data
        assert right.base is standard.data

    def test_amp33(self, data):
        standard = StandardView(data)
        amp33 = standard.amp33

        # Check that the amp33 is the correct shape and type
        assert amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.CHAN_WIDTH)
        assert amp33.dtype == np.float32

        # Check that the amp33 is the correct part of the data
        assert (amp33 == data[:, :, -Const.CHAN_WIDTH :]).all()

        # Check that the amp33 is a view of the data
        assert amp33.base is data
        assert amp33.base is standard.data

    def test_data_view_relations(self, standard):
        """
        Check the relationships among the data views.
        """
        # Check that the detector + amp33 is the same as the data
        assert (np.dstack([standard.detector, standard.amp33]) == standard.data).all()

        # Check that the left is the left side of the detector
        assert (standard.left == standard.detector[:, :, : Const.REF]).all()

        # Check that the right is the right side of the detector
        assert (standard.right == standard.detector[:, :, -Const.REF :]).all()

    def test_channels(self, data, offset):
        """
        Check aligning the data into channels and padding the columns.
        """
        non_view_data = data.copy()
        standard = StandardView(data, offset=offset)
        channel_view = standard.channels

        # Check the data is unchanged by the property
        assert (standard.data == non_view_data).all()

        # Check the channel_view is a ChannelView object
        assert isinstance(channel_view, ChannelView)
        assert channel_view.offset is standard.offset

        channels = channel_view.data

        # check that the data is not a view
        assert channels.base is None

        # Check that the channels are the correct shape and type
        assert channels.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Const.CHAN_WIDTH + Const.PAD,
        )
        assert channels.dtype == np.float32

        # Check the channels
        for idx in range(0, Dims.N_CHAN):
            truth = data[:, :, idx * Const.CHAN_WIDTH : (idx + 1) * Const.CHAN_WIDTH]
            if idx % 2 == 0:
                # Even channels are not reversed
                assert (channels[idx, :, :, : Const.CHAN_WIDTH] == truth).all()
            else:
                # Odd channels are reversed
                assert (
                    channels[idx, :, :, : Const.CHAN_WIDTH] == truth[:, :, ::-1]
                ).all()

        # Check that amp33 is the last channel
        assert (standard.amp33 == channels[-1, :, :, : Const.CHAN_WIDTH]).all()

        # Check that the first part of the first channel is left
        assert (standard.left == channels[0, :, :, : Const.REF]).all()

        # Check that the last part of the last detector channel is right
        #     This is an odd indexed channel, so it is reversed, this also means
        #     the first 4 columns in the channel form are the ones corresponding
        #     to the right side of the detector (not the last 4 columns) as one
        #     might expect.
        assert (standard.right[:, :, ::-1] == channels[-2, :, :, : Const.REF]).all()

    def test_channels_regression(self, data):
        standard = StandardView(data)
        regression = data.copy()

        # Run the regression utility
        regression = reference_utils.aligned_channels(regression)

        # Check that the regression matches the channels property
        assert (standard.channels.data == regression).all()

    def test_remove_offset(self, data):
        non_view_data = data.copy()
        standard = StandardView(data)
        assert standard.offset is None

        new = standard.remove_offset()

        # Check that new object returned is the same as the original
        assert new is standard

        # Check that the data is still a view of the original
        assert new.data is data

        # Check that the offset is now set
        assert new.offset is not None
        assert new.offset.shape == (Dims.N_ROWS, Dims.N_COLS)
        assert (new.offset != 0).all()

        # Check that the data has been updated and that the offset has been removed
        assert (new.data != non_view_data).any()
        assert (new.data == (non_view_data - new.offset)).all()

        # add offset back (error accumulates due to floating point arithmetic)
        reset = StandardView(new.data + new.offset)
        assert_allclose(reset.data, non_view_data, atol=1e-5)

    def test_remove_offset_regression(self, data):
        standard = StandardView(data)

        # Copy the data, because the regression utility modifies the data in-place
        regression = data.copy()
        _, b = reference_utils.remove_linear_trends(regression, True)

        # Run the internal utility
        new = standard.remove_offset()

        # Check that the regression matches the new object
        assert (new.offset == b).all()
        assert (new.data == regression).all()

    def test_apply_offset(self, data, offset):
        standard = StandardView(data.copy(), offset)

        # Apply the offset
        new = standard.apply_offset()

        # Check the new object is the same as the original
        assert new is standard

        # Check apply offset has modified new's state
        assert new.offset is None
        assert (new.data != data).any()

        # Check that the data has been correctly modified
        for new_frame, frame in zip(new.data, data):
            assert (new_frame == (frame + offset)).all()

        # Check if no offset is passed
        standard = StandardView(data.copy())
        new = standard.apply_offset()
        assert new is standard
        assert (new.data == data).all()

    def test_apply_offset_regression(self, data, offset):
        standard = StandardView(data.copy(), offset)

        # Run the reference utility
        regression = data[:, :, : Const.N_COLUMNS].copy()
        reference_utils.restore_offsets(Dims.N_FRAMES, regression, offset)

        # Run the internal utility
        new = standard.apply_offset()

        # Check regression
        assert (new.data[:, :, : Const.N_COLUMNS] == regression).all()

    def test_offset_roundtrip(self, data):
        standard = StandardView(data.copy())

        no_offset = standard.remove_offset()
        assert no_offset.offset is not None
        assert (no_offset.data != data).any()

        reset = no_offset.apply_offset()
        assert reset.offset is None

        # Check round trip has not changed the data.
        #   The differences are due to floating point arithmetic.
        assert_allclose(reset.data, data, atol=1e-5)


class TestChannelView:
    def test_detector(self, standard):
        channels = standard.channels
        detector = channels.detector

        # Check that the detector is the correct shape and type
        assert detector.shape == (
            Const.N_DETECT_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.PADDED_WIDTH,
        )
        assert detector.dtype == np.float32

        # Check that the detector is the correct part of the data
        assert (detector == channels.data[: Const.N_DETECT_CHAN, :, :, :]).all()

        # Check that the detector is a view of the data
        assert detector.base is channels.data

    def test_left(self, standard):
        channels = standard.channels
        left = channels.left

        # Check that the left is the correct shape and type
        assert left.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.PADDED_WIDTH)
        assert left.dtype == np.float32

        # Check the padding is correct
        assert (left[:, :, Const.REF :] == 0).all()

        # Check that the left is the correct part of the data
        assert (left[:, :, : Const.REF] == channels.data[0, :, :, : Const.REF]).all()
        assert (
            left[:, :, : Const.REF] == channels.detector[0, :, :, : Const.REF]
        ).all()
        assert (left[:, :, : Const.REF] == standard.left).all()

        # Check that the left is not a view
        assert left.base is None

    def test_right(self, standard):
        channels = standard.channels
        right = channels.right

        # Check that the right is the correct shape and type
        assert right.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.PADDED_WIDTH)
        assert right.dtype == np.float32

        # Check the padding is correct
        assert (right[:, :, Const.REF :] == 0).all()

        # Check that the right is the correct part of the data
        assert (right[:, :, : Const.REF] == channels.data[-2, :, :, : Const.REF]).all()
        assert (
            right[:, :, : Const.REF] == channels.detector[-1, :, :, : Const.REF]
        ).all()

        # Due to the channel flip the right will need to be reversed to match the
        #   standard.right
        assert (right[:, :, : Const.REF] == standard.right[:, :, ::-1]).all()

        # Check that the right is not a view
        assert right.base is None

    def test_amp33(self, standard):
        channels = standard.channels
        amp33 = channels.amp33

        # Check that the right is the correct shape and type
        assert amp33.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.PADDED_WIDTH)
        assert amp33.dtype == np.float32

        # Check that the right is the correct part of the data
        assert (amp33 == channels.data[-1, :, :, :]).all()
        assert (amp33[:, :, : Const.CHAN_WIDTH] == standard.amp33).all()

        # Check that the amp33 is a view
        assert amp33.base is channels.data

    def test_standard(self, standard, offset):
        channels = standard.channels
        channels.offset = offset
        non_view_data = channels.data.copy()

        standard_view = channels.standard

        # Check the channels.data is unchanged by the property
        assert (channels.data == non_view_data).all()

        # Check the standard_view is a StandardView object
        assert isinstance(standard_view, StandardView)
        assert standard_view.offset is channels.offset

        data = standard_view.data

        # Check the data is not a view of the channels.data
        assert data.base is not channels.data

        # Check that the data is the correct shape and type
        assert data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)
        assert data.dtype == np.float32

        # Check that data matches the standard data (round trip)
        assert (data == standard.data).all()

    def test_remove_trends(self, standard):
        channels = standard.channels
        non_view_data = channels.data.copy()

        new = channels.remove_trends()

        # Check that new object returned is the original
        assert new is channels

        # Check the data is not a view
        assert new.data.base is None

        # Check the method has updated the data
        assert (new.data != non_view_data).any()

    def test_remove_trends_regression(self, standard):
        channels = standard.channels

        # Copy the data, because the regression utility modifies the data in-place
        regression = channels.data.copy()
        reference_utils.remove_linear_trends_per_frame(regression, False, True)

        # Run the internal utility
        new = channels.remove_trends()

        # Check that the regression matches the new object
        assert (new.data == regression).all()

    def test_cosine_interpolate(self, standard):
        channels = standard.channels
        non_view_data = channels.data.copy()

        new = channels.cosine_interpolate()

        # Check that new object returned is the original
        assert new is channels

        # Check the data is not a view
        assert new.data.base is None

        # Check the method has updated only the amp33 data
        assert (new.amp33 != non_view_data[-1, :, :, :]).any()
        assert (new.detector == non_view_data[:-1, :, :, :]).all()

    def test_cosine_interpolate_regression(self, channels):
        # Copy the data, because the regression utility modifies the data in-place
        regression = channels.data.copy()
        reference_utils.cos_interp_reference(regression, regression.shape[1])

        # Demonstrate the regression only ran on the amp33 data
        assert (regression[-1, :, :, :] != channels.data[-1, :, :, :]).any()
        assert (regression[:-1, :, :, :] == channels.data[:-1, :, :, :]).all()

        # Run the internal utility
        new = channels.cosine_interpolate()

        # Check that the regression matches the new object
        assert (new.amp33 == regression[-1, :, :, :]).all()
        assert (new.data == regression).all()

    def test_fft_interpolate(self, channels):
        non_view_data = channels.data.copy()

        new = channels.fft_interpolate()

        # Check that the new object returned is the original
        assert new is channels

        # Check the data is not a view
        assert new.data.base is None

        # Check the method has updated only the amp33 data
        assert (new.amp33 != non_view_data[-1, :, :, :]).any()
        assert (new.detector == non_view_data[:-1, :, :, :]).all()

    def test_fft_interpolate_regression(self, channels):
        """
        Run fft interpolation regression test
        """
        non_view_data = channels.data.copy()

        # Copy the data, because the regression utility modifies the data in-place
        regression = channels.data.copy()
        amp33_regression = reference_utils.fft_interp_amp33(
            regression, regression.shape[1]
        )

        # Demonstrate that the sub-array (amp33_regression) does have changes made
        assert (amp33_regression[:, :, :] != non_view_data[-1, :, :, :]).any()

        # Demonstrate that the actual regression "data" remains unchanged instead
        # the claimed in-place change. (this would fail if in-place change was made)
        assert (regression == non_view_data).all()

        # Run the internal utility
        new = channels.fft_interpolate()

        assert (new.amp33 == amp33_regression).all()

    def test_reference_fft(self, channels):
        reference = channels.reference_fft

        # Check the new object is a ReferenceFFT object
        assert isinstance(reference, ReferenceFFT)

        # Check the FFT all have the correct shape
        shape = (Dims.N_FRAMES, (Dims.N_ROWS * Dims.PADDED_WIDTH) // 2 + 1)
        assert reference.left.shape == shape
        assert reference.right.shape == shape
        assert reference.amp33.shape == shape

        # Check the FFT are all complex dtype
        assert reference.left.dtype == np.complex64
        assert reference.right.dtype == np.complex64
        assert reference.amp33.dtype == np.complex64

        # Check the FFT are not views
        assert reference.left.base is None
        assert reference.right.base is None
        assert reference.amp33.base is None

    def test_reference_fft_regression(self, channels):
        regression = channels.data.copy()
        reference = channels.reference_fft

        # Sow that the reference property does not alter the data
        assert (regression == channels.data).all()

        # run the reference code
        left, right, amp33 = reference_utils.forward_fft(Dims.N_FRAMES, regression)

        # Copy to make sure regression utility does not modify the data in-place
        assert (regression == channels.data).all()

        # Check that the regression matches the new object
        assert (reference.left == left).all()
        assert (reference.right == right).all()
        assert (reference.amp33 == amp33).all()

    def test_correction(self, channels, coeffs):
        correction = channels.correction(coeffs)

        # Check size and dtype
        assert correction.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.PADDED_WIDTH,
        )
        assert correction.dtype == np.float32

    def test_correction_regression(self, channels, coeffs):
        # Run the reference code
        regression = reference_utils.correction(
            Dims.N_FRAMES, channels.data.copy(), coeffs.alpha, coeffs.gamma, coeffs.zeta
        )

        # Check size and dtype
        assert regression.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Dims.PADDED_WIDTH,
        )
        assert regression.dtype == np.float64

        # Run the internal utility
        correction = channels.correction(coeffs)

        assert (correction == regression.astype(np.float32)).all()

    def test_apply_correction(self, channels, coeffs):
        non_view_standard_data = channels.standard.data.copy()
        non_view_data = channels.data.copy()

        amp33 = channels.amp33.copy()

        # Compute the applied correction
        standard = channels.apply_correction(coeffs)

        # Check the correction is the same as the standard
        assert isinstance(standard, StandardView)

        # Check the channels is updated in-place
        assert (channels.data != non_view_data).any()

        # Check the standard is different from the original
        assert (standard.data != non_view_standard_data).any()

        # Check the amp33 channel has not changed
        assert (channels.amp33 == amp33).all()

    def test_apply_correction_channels_regression(self, channels, coeffs):
        # Get regression
        regression = reference_utils.apply_correction_to_channels(
            Dims.N_FRAMES, channels.data.copy(), coeffs.alpha, coeffs.gamma, coeffs.zeta
        )

        # Check size and dtype
        assert regression.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS,
            Const.CHAN_WIDTH,
        )
        assert regression.dtype == np.float32

        # Run the internal utility (looking at the inplace changes)
        channels.apply_correction(coeffs)

        # Check regression
        assert (channels.data[:, :, :, : Const.CHAN_WIDTH] == regression).all()

    def test_apply_correction_standard_regression(self, channels, coeffs):
        regression = reference_utils.apply_correction_to_data(
            Dims.N_FRAMES, channels.data.copy(), coeffs.alpha, coeffs.gamma, coeffs.zeta
        )

        # Check size and dtype
        assert regression.shape == (Dims.N_FRAMES, Dims.N_ROWS, Const.N_COLUMNS)
        assert regression.dtype == np.float32

        # Run the internal utility
        standard = channels.apply_correction(coeffs)

        # Check regression
        assert (standard.detector == regression).all()


class TestReferenceFFT:
    def test_correction(self, channels, coeffs):
        reference = channels.reference_fft

        left = reference.left.copy()
        right = reference.right.copy()
        amp33 = reference.amp33.copy()

        correction = reference.correction(coeffs)

        # Check size and dtype
        assert correction.shape == (
            Dims.N_CHAN,
            Dims.N_FRAMES,
            Dims.N_ROWS * Dims.PADDED_WIDTH,
        )
        assert correction.dtype == np.float64  # no dtype casting yet

        # Check the last channel is 0
        assert (correction[-1, :, :] == 0).all()

        # Check the FFT have not been modified
        assert (left == reference.left).all()
        assert (right == reference.right).all()
        assert (amp33 == reference.amp33).all()

    def test_correction_regression(self, channels, coeffs):
        reference = channels.reference_fft

        # Run reference code
        regression = reference_utils.compute_correction(
            Dims.N_ROWS,
            coeffs.alpha,
            coeffs.gamma,
            coeffs.zeta,
            Dims.N_FRAMES,
            reference.left,
            reference.right,
            reference.amp33,
        )

        # Run implementation
        correction = reference.correction(coeffs)

        assert correction.dtype == regression.dtype
        assert (correction == regression).all(), correction.dtype


class TestCoefficients:
    def test_from_ref(self, ref_pix_ref):
        coeffs = Coefficients.from_ref(ref_pix_ref)

        # Check the object is a Coefficients object
        assert isinstance(coeffs, Coefficients)

        # Check the coefficients are the same
        assert (coeffs.alpha == ref_pix_ref.alpha).all()
        assert (coeffs.gamma == ref_pix_ref.gamma).all()
        assert (coeffs.zeta == ref_pix_ref.zeta).all()
