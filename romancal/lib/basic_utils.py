"""General utility objects"""

import numpy as np


def parse_visitID(visit_id):
    """Utility to parse the visit_id into its components

    Input:
       visit_id as a string

    Output:
      program number
      execution plan number
      pass number
      segment number
      observation number
      visit number
    """

    visit_id_parts = dict(
        [
            ("Program", visit_id[0:5]),
            ("Execution", visit_id[5:7]),
            ("Pass", visit_id[7:10]),
            ("Segment", visit_id[10:13]),
            ("Observation", visit_id[13:16]),
            ("Visit", visit_id[16:20]),
        ]
    )

    return visit_id_parts


def frame_read_times(frame_time, sca, frame_number=0):
    """
    Compute the pixel read times for a single frame.

    This is a placeholder function that assumes a uniform read
    across each channel within the frame time.  A more careful
    treatment will need to account for time spent reading out the
    guide window.

    Data shape for the frame is assumed to be 4096 x 4096, with
    32 channels along the columns.

    Parameters
    ----------
    frame_time : float
        The frame time for the exposure, in seconds.
    sca : int
        The WFI detector number (1-18).
    frame_number : float, optional
        The frame number.  Default of zero means that pixels start
        reading out at t=0.

    Returns
    -------
    read_times : `~numpy.ndarray`
        A 4096 x 4096 array containing the read time in seconds for each pixel.
    """
    nchannel = 32
    nrow = 4096
    ncol = 128
    one_channel_read_time = np.linspace(
        0, frame_time, nrow * ncol, endpoint=False
    ).reshape(nrow, ncol)

    # WFI channels alternate readout direction in the +x and -x directions
    # we implement this by flipping the x direction of every other channel
    two_channel_read_times = np.concatenate(
        [one_channel_read_time, one_channel_read_time[:, ::-1]], axis=1
    )
    read_times = np.tile(two_channel_read_times, (1, nchannel // 2))

    # Apply science -> detector flipping for read order.
    # Detectors with SCA % 3 == 0 flip columns; others flip rows.
    if sca % 3 == 0:
        read_times = read_times[:, ::-1]
    else:
        read_times = read_times[::-1, :]

    read_times += frame_number * frame_time

    return read_times


def compute_var_rnoise(model):
    """Compute read noise variance from model data.

    If var_rnoise exists in the model, return it directly.
    Otherwise, compute it as err^2 - sum(other variance terms).

    This function supports the optional storage of var_rnoise in L2 files.
    When var_rnoise is not present, it can be reconstructed from the total
    error and other variance components using:
        var_rnoise = err^2 - var_poisson - var_flat - var_dark - ...

    Parameters
    ----------
    model : ImageModel
        Roman WFI ImageModel containing error and variance arrays.

    Returns
    -------
    var_rnoise : np.ndarray
        Read noise variance array.

    Notes
    -----
    The total error follows the relation:
        err^2 = var_rnoise + var_poisson + var_flat + var_dark + ...

    Therefore:
        var_rnoise = err^2 - var_poisson - var_flat - var_dark - ...
    """
    # If var_rnoise exists in the model, return it
    if hasattr(model, "var_rnoise"):
        return model.var_rnoise

    # Otherwise, compute from err^2 minus other variance terms
    var_rnoise = model.err.astype(np.float32) ** 2

    # Subtract other variance components
    variance_arrays = ["var_poisson", "var_flat", "var_dark"]
    for var_name in variance_arrays:
        if hasattr(model, var_name):
            var_rnoise -= getattr(model, var_name)

    return var_rnoise
