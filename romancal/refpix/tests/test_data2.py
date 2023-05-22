from enum import IntEnum

import numpy as np
import pytest
from numpy.testing import assert_allclose

from romancal.refpix.data2 import ChannelView, Const, StandardView

from . import reference_utils

RNG = np.random.default_rng(42)


class Dims(IntEnum):
    """
    Assumed data dimensions for dat (aid clarity of testing).
        - N_FRAMES: the number of frames used in the test data (arbitrary) choice
        - N_ROWS: the number of rows of the test data
        - N_COLS: the number of columns of the test data
    """

    N_FRAMES = 8
    N_ROWS = Const.N_COLUMNS
    N_COLS = Const.N_COLUMNS + Const.CHAN_WIDTH
    N_CHAN = N_COLS // Const.CHAN_WIDTH
    PADDED_WIDTH = Const.CHAN_WIDTH + Const.PAD


@pytest.fixture(scope="module")
def data():
    return RNG.uniform(1, 100, size=(Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)).astype(
        np.float32
    )


@pytest.fixture(scope="module")
def standard(data):
    return StandardView(data)


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
    assert Const.N_COLUMNS == 4096 == Dims.N_ROWS
    assert Dims.N_COLS == 4224

    # Check that there are the correct total number of channels
    assert Dims.N_CHAN == 33 == Const.N_DETECT_CHAN + 1

    # Check the padded channel width
    assert Dims.PADDED_WIDTH == 140 == Const.CHAN_WIDTH + Const.PAD


def test_data(data):
    """
    Sanity that the generated data is the correct shape and type.
    """
    assert data.shape == (Dims.N_FRAMES, Dims.N_ROWS, Dims.N_COLS)
    assert data.dtype == np.float32


class TestStandardView:
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

    def test_channels(self, data):
        """
        Check aligning the data into channels and padding the columns.
        """
        non_view_data = data.copy()
        standard = StandardView(
            data, offset=RNG.uniform(0, 100, size=(Dims.N_ROWS, Dims.N_COLS))
        )
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

    def test_standard(self, standard):
        channels = standard.channels
        channels.offset = RNG.uniform(0, 100, size=(Dims.N_ROWS, Dims.N_COLS))
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

    def test_cosine_interpolate_regression(self, standard):
        channels = standard.channels

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

    def test_fft_interpolate(self, standard):
        channels = standard.channels
        non_view_data = channels.data.copy()

        new = channels.fft_interpolate()

        # Check that the new object returned is the original
        assert new is channels

        # Check the data is not a view
        assert new.data.base is None

        # Check the method has updated only the amp33 data
        assert (new.amp33 != non_view_data[-1, :, :, :]).any()
        assert (new.detector == non_view_data[:-1, :, :, :]).all()

    def test_fft_interpolate_regression(self, standard):
        """
        Run fft interpolation regression test

        NOTE:
            The reference code assumes the data will be changed in-place for all
            its major operations.However, the reference code violates this assumption
            for the Amp33 FFT interpolation step. It does make an in-place change to
            a sub-array, `dataReferenceChan_FramesFlat`. However
            `dataReferenceChan_FramesFlat` loses its view on the original data array
            because it unnecessarily recasts the dtype, which is a non-view compatible
            change.

            romancal does not make this mistake because it does not recast the dtype.

            The reference utils are modified to output the data array which is modified
            as a functional return value. This is so that the reference code's
            computation can be tests against the romancal code's computation.

            After this step is computed, the reference code an romancal code will
            diverge in values, because the changes made by this will propagate in
            romancal but the do not in the reference code.
        """

        channels = standard.channels
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
