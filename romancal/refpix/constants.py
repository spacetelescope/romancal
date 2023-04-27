"""
Constant values for reference pixel correction.
"""

from enum import IntEnum


class RefPix(IntEnum):
    """Reference pixel types."""

    LEFT_CHANNEL = 0
    RIGHT_CHANNEL = 31
    REF_CHANNEL = 32
    CHANNEL_WIDTH = 4


# Number of rows in the detector
NUM_ROWS = 4096

# Fixed channel values for Roman
NUM_CHANNELS = 32

NUM_COLUMNS_PER_CHANNEL = 128
NUM_COLUMNS = NUM_CHANNELS * NUM_COLUMNS_PER_CHANNEL

# The effective number of pixels sampled during the pause at the end of each read
# this preserves the phase of temporally periodic signals
END_PAD = 12
NUM_COLUMNS_PER_CHANNEL_PAD = NUM_COLUMNS_PER_CHANNEL + END_PAD

# Data size per channel
SIZE_PER_CHANNEL = NUM_COLUMNS_PER_CHANNEL_PAD * NUM_ROWS

# Reference pixel channels
LEFT_CHANNEL = 0
RIGHT_CHANNEL = 31
REF_CHANNEL = 32

# The left and right reference pixel streams are 4 pixels wide
CHANNEL_WIDTH = 4

# The reordering of the channel representation,
#    Original: (frame, row, channel, channel_column)
CHANNEL_REORDER = (2, 0, 1, 3)  # Reorders to: (channel, frame, row, channel_column)
CHANNEL_RESET = (1, 2, 0, 3)  # Undo the reordering

# Stream normalization factor for left and right reference pixel streams
STREAM_NORM = float(NUM_COLUMNS_PER_CHANNEL_PAD) / float(CHANNEL_WIDTH)
