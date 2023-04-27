# import numpy as np
# from scipy import fft

# from . import constants as const


# def remove_offset(data, amp33):


# def split_channels(data):
#     """
#     Split the input data into channels for each amplifier.
#     """

#     num_frames = data.shape[0]

#     # Split data into channels [frame, row, channel, channel_column]
#     main_data = data.reshape((num_frames, const.NUM_ROWS, const.NUM_CHANNELS,
#                               const.CHANNEL_WIDTH))

#     # Reorder the channels to align data so computation happens on contiguous memory
#     #    [channel, frame, row, channel_column]
#     main_data = data.transpose(const.CHANNEL_REORDER)

#     # Reverse every other channel column so columns are in reading order
#     for channel in range(1, const.NUM_CHANNELS, 2):
#         data[channel, :, :, :] = main_data[channel, :, :, ::-1]

#     # pad channels with zeros to account for the pause at the end of each read
#     return np.pad(main_data, ((0, 0), (0, 0), (0, 0), (0, const.END_PAD)),
#                   mode='constant', constant_values=0)


# def extract_reference_streams(data):
#     """Extract the reference columns from the input data.

#     Parameters
#     ----------
#     data : numpy.ndarray
#         The input data.

#     Returns
#     -------
#     left, right, ref : numpy.ndarray
#         The left, right, and reference streams.
#     """

#     num_frames = data.shape[0]

#     def _extract(channel):
#         return np.copy(data[channel, :, :, :])

#     def _clear(stream):
#         stream[:, :, const.CHANNEL_WIDTH:] = 0


#     left = _extract(const.LEFT_CHANNEL)
#     _clear(left)

#     right = _extract(const.RIGHT_CHANNEL)
#     _clear(right)

#     ref = np.copy(data[const.REF_CHANNEL, :, :, :])

#     def _reshape(stream):
#         return stream.reshape((num_frames, const.SIZE_PER_CHANNEL))

#     return _reshape(left), _reshape(right), _reshape(ref)


# def fft_reference_streams(left, right, ref):
#     def _fft(stream):
#         return fft.rfft(stream / stream[0].size)

#     return const.STREAM_NORM * _fft(left), const.STREAM_NORM * _fft(right), _fft(ref)
