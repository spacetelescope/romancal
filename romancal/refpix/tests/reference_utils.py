import threading

import numpy as np

NUM_ROWS = 4096
NUM_OUTPUT_CHANS = 33
NUM_COLS_PER_OUTPUT_CHAN = 128
NUM_COLS = NUM_OUTPUT_CHANS * NUM_COLS_PER_OUTPUT_CHAN
END_OF_ROW_PIXEL_PAD = 12
NUM_COLS_PER_OUTPUT_CHAN_WITH_PAD = NUM_COLS_PER_OUTPUT_CHAN + END_OF_ROW_PIXEL_PAD
ALL_CHAN_RANGE = range(NUM_OUTPUT_CHANS)
REFERENCE_ROWS = [0, 1, 2, 3, 4092, 4093, 4094, 4095]


def remove_linear_trends(data_FramesRowsCols, subtractOffsetOnly):
    """
    Remove linear trends (gains/slope, biases/offsets) per pixel
    :param data_FramesRowsCols: image data [numFrames][numRows][numCols].
        The data is modified IN PLACE
    :param subtractOffsetOnly: if True, only the liner model offset is removed.
        If False, the offset and slope are removed
    :return per-pixel image array of modeled m (slope) and b (offset)
    """

    numFrames = len(data_FramesRowsCols)
    numRows = len(data_FramesRowsCols[0])
    numCols = len(data_FramesRowsCols[0][0])

    # for mat multiply
    data_FramesRowsCols = data_FramesRowsCols.reshape((numFrames, numCols * numRows))

    frameRamp = np.arange(numFrames, dtype=data_FramesRowsCols.dtype)

    # subtract mean(frameRampMinusMean) to minimize correlation of slope and offset
    frameRampMinusMean = frameRamp - np.mean(frameRamp)

    sxy = np.matmul(
        frameRampMinusMean, data_FramesRowsCols, dtype=data_FramesRowsCols.dtype
    )
    sx = np.sum(frameRampMinusMean)
    sxx = np.sum(frameRampMinusMean**2)
    sy = np.sum(data_FramesRowsCols, axis=0)

    m = (numFrames * sxy - sx * sy) / (numFrames * sxx - sx**2)
    b = (sy * sxx - sxy * sx) / (numFrames * sxx - sx**2)

    for frame in range(numFrames):
        if not subtractOffsetOnly:
            data_FramesRowsCols[frame, :] -= frameRampMinusMean[frame] * m + b
        else:
            data_FramesRowsCols[frame, :] -= b  # subtract offset

    data_FramesRowsCols = data_FramesRowsCols.reshape(numFrames, numRows, numCols)
    return m.reshape(numRows, numCols), b.reshape(numRows, numCols)


def aligned_channels(data):
    numFrames = data.shape[0]

    data_chansFramesRowsPhyscols = np.transpose(
        data.reshape((numFrames, NUM_ROWS, NUM_OUTPUT_CHANS, NUM_COLS_PER_OUTPUT_CHAN)),
        (2, 0, 1, 3),
    )

    for chanNum in range(1, NUM_OUTPUT_CHANS, 2):
        data_chansFramesRowsPhyscols[chanNum, :, :, :] = data_chansFramesRowsPhyscols[
            chanNum, :, :, ::-1
        ]

    return np.pad(
        data_chansFramesRowsPhyscols,
        ((0, 0), (0, 0), (0, 0), (0, END_OF_ROW_PIXEL_PAD)),
    )


def remove_linear_trends_per_frame(
    data_chansFramesRowsChancols, subtractOffsetOnly, multiThread=True
):
    """
    Entry point for Fitting and removal of slopes per frame to remove issues
        at frame boundaries.
    :param data_chansFramesRowsChancols:
        [numOutChannels(33)][numFrames][numRows][numColsWithPadding].
        The data is modified IN PLACE
    :param subtractOffsetOnly: if True, only the liner model offset is removed.
        If False, the offset and slope are removed
    """

    numFrames = len(data_chansFramesRowsChancols[0])
    numRows = len(data_chansFramesRowsChancols[0][0])
    timeDomainRange = np.arange(
        NUM_COLS_PER_OUTPUT_CHAN_WITH_PAD * numRows,
        dtype=data_chansFramesRowsChancols.dtype,
    ).reshape(numRows, NUM_COLS_PER_OUTPUT_CHAN_WITH_PAD)
    timeDomainRange = timeDomainRange - timeDomainRange.mean()

    exec_channel_func_threads(
        ALL_CHAN_RANGE,
        remove_linear_trends_per_frame_chan_func,
        (numFrames, timeDomainRange, data_chansFramesRowsChancols, subtractOffsetOnly),
        multiThread=multiThread,
    )


def remove_linear_trends_per_frame_chan_func(
    chanNumber,
    numFrames,
    timeDomainRange,
    data_chansFramesRowsChancols,
    subtractOffsetOnly,
):
    """
    Per-channel function for fitting and removal of slopes per frame to remove
        issues at frame boundaries - single channel only
    :param chanNumber:
    :param numFrames:
    :param timeDomainRange:
    :param data_chansFramesRowsChancols:
        [numOutChannels(33)][numFrames][numRows][numColsWithPadding].
        The data is modified IN PLACE
    :param subtractOffsetOnly: if True, only the liner model offset is removed.
        If False, the offset and slope are removed
    """

    ab_3 = np.zeros(
        (NUM_OUTPUT_CHANS, numFrames, 2), dtype=np.float64
    )  # for top+bottom ref pixel rows
    tdRef = timeDomainRange[REFERENCE_ROWS, :]

    for frame in range(numFrames):
        xVals = tdRef[
            (data_chansFramesRowsChancols[chanNumber, frame, REFERENCE_ROWS, :] != 0)
        ]  # where data[chanNumber,frame,row4plus4,:] != 0)
        data2 = data_chansFramesRowsChancols[chanNumber, frame, REFERENCE_ROWS, :]
        yVals = data2[(data2 != 0)]
        if xVals.size < 1 or yVals.size < 1:
            continue

        ab_3[chanNumber, frame, :] = np.polyfit(xVals, yVals, 1)

        dataNotZero = data_chansFramesRowsChancols[chanNumber, frame, :, :] != 0

        if not subtractOffsetOnly:
            data2 = (
                timeDomainRange * ab_3[chanNumber, frame, 0]
                + ab_3[chanNumber, frame, 1]
            )
        else:
            data2 = ab_3[chanNumber, frame, 1]

        data_chansFramesRowsChancols[chanNumber, frame, :, :] -= data2 * dataNotZero


def exec_channel_func_threads(chanIndexRange, targetFunc, funcArgs, multiThread=False):
    """
    Execute a function over a range of channels.  If multiThread is True, all executions
    will be executed on individual threads.

    :param chanIndexRange: E.g., 'range(NUM_OUTPUT_CHANS')
    :param targetFunc: function to call, first arg must be channel number
    :param funcArgs: function argument tuple NOT INCLUDING channel number e.g.
        (sig_data, goodPixelsOneIsMask_RowCol)
    :param multiThread: if True allocate a separate thread for each channel;
        otherwise the channels are processed sequentially in a single thread

    """

    if multiThread:
        threadList = []
        for c in chanIndexRange:
            funcArgsWithChan = (c,) + funcArgs
            threadList.append(
                threading.Thread(target=targetFunc, args=funcArgsWithChan)
            )

        for t in threadList:
            t.start()

        for t in threadList:
            t.join()
    else:
        for c in chanIndexRange:
            funcArgsWithChan = (c,) + funcArgs
            targetFunc(*funcArgsWithChan)
