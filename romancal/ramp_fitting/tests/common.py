import numpy as np
from astropy.time import Time
from roman_datamodels.datamodels import (
    DarkRefModel,
    GainRefModel,
    RampModel,
    ReadnoiseRefModel,
)

# Read Time in seconds
#   For Roman, the read time of the detectors is a fixed value and is currently
#   backed into code. Will need to refactor to consider the more general case.
#   Used to deconstruct the MultiAccum tables into integration times.
ROMAN_READ_TIME = 3.04

rng = np.random.default_rng(42)

# Basic resultant
#
# The read pattern is `[[1], [2], [3], [4]]`
# The total expected counts is 7.
# The resultants were generated with
# `romanisim.l1.apportion_counts_to_resultants(counts, read_pattern)`.
SIMPLE_RESULTANTS = np.array(
    [
        [[2.0, 2.0], [5.0, 1.0]],
        [[4.0, 5.0], [6.0, 2.0]],
        [[5.0, 6.0], [7.0, 6.0]],
        [[7.0, 7.0], [7.0, 7.0]],
    ],
    dtype=np.float32,
)


# #########
# Utilities
# #########
def make_data(resultants, ingain, rnoise, randomize):
    """Create test input data

    Parameters
    ----------
    resultants : numpy.array.shape(3, xdim, ydim)
        The resultant array.

    ingain : int
        Gain to apply

    rnoise : float
        Noise to apply

    randomize : bool, expected)
        Randomize the gain and read noise across pixels.

    Returns
    -------
    image, gain, readnoise : ImageModel, GainRefModel, ReadnoiseRefModel.array
        Input image and related references
    """
    ramp_model = model_from_resultants(resultants)
    gain_model, readnoise_model, dark_model = generate_wfi_reffiles(
        ramp_model.shape[1:], ingain=ingain, rnoise=rnoise, randomize=randomize
    )

    return ramp_model, gain_model, readnoise_model, dark_model


def model_from_resultants(resultants, read_pattern=None):
    """Create a RampModel from resultants

    Parameters
    ----------
    resultants : numpy.array.shape(reads, xdim, ydim)
        The resultants to fit.

    read_pattern : [[int[,...]][,...]]
        The read pattern used to produce the resultants.
        If None, presume a basic read pattern
    """
    if read_pattern is None:
        read_pattern = [[idx + 1] for idx in range(resultants.shape[0])]

    # Full WFI image has reference pixels all around. Add those on.
    nrefpixs = 4
    full_wfi = np.ones(
        (
            resultants.shape[0],
            resultants.shape[1] + (nrefpixs * 2),
            resultants.shape[2] + (nrefpixs * 2),
        ),
        dtype=np.float32,
    )
    full_wfi[:, nrefpixs:-nrefpixs, nrefpixs:-nrefpixs] = resultants
    shape = full_wfi.shape

    pixdq = np.zeros(shape=shape[1:], dtype=np.uint32)
    err = np.zeros(shape=shape, dtype=np.float32)
    gdq = np.zeros(shape=shape, dtype=np.uint8)

    ramp_model = RampModel.create_fake_data(shape=shape)
    ramp_model.meta.cal_step = {}
    for step_name in ramp_model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        ramp_model.meta.cal_step[step_name] = "INCOMPLETE"
    ramp_model.data = full_wfi
    ramp_model.pixeldq = pixdq
    ramp_model.groupdq = gdq
    ramp_model.err = err

    ramp_model.meta.exposure.frame_time = ROMAN_READ_TIME
    ramp_model.meta.exposure.nresultants = shape[0]

    ramp_model.meta.exposure.read_pattern = read_pattern

    return ramp_model


def generate_wfi_reffiles(
    shape, ingain=6, rnoise=0.01, darkcurrent=0.01, randomize=True
):
    """Create GainRefModel and ReadnoiseRefModel

    Parameters
    ----------
    shape : tuple
        Shape of the arrays

    ingain : float
        Maximum gain.

    rnoise : flota
        Maximum noise

    darkcurrent : float
        Dark current scale

    randomize : bool
        Randomize the gain and read noise data.
    """
    # Create temporary gain reference file
    gain_ref_model = GainRefModel.create_fake_data(shape=shape)

    gain_ref_model["meta"]["instrument"]["detector"] = "WFI01"
    gain_ref_model["meta"]["instrument"]["name"] = "WFI"
    gain_ref_model["meta"]["reftype"] = "GAIN"
    gain_ref_model["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    if randomize:
        gain_ref_model["data"] = (rng.random(shape) * 0.5).astype(np.float32) * ingain
    else:
        gain_ref_model["data"] = np.ones(shape).astype(np.float32) * ingain
    gain_ref_model["dq"] = np.zeros(shape, dtype=np.uint16)
    gain_ref_model["err"] = (rng.random(shape) * 0.05).astype(np.float32)

    # Create temporary readnoise reference file
    rn_ref_model = ReadnoiseRefModel.create_fake_data(shape=shape)
    rn_ref_model["meta"]["instrument"]["detector"] = "WFI01"
    rn_ref_model["meta"]["instrument"]["name"] = "WFI"
    rn_ref_model["meta"]["reftype"] = "READNOISE"
    rn_ref_model["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    rn_ref_model["meta"]["exposure"]["type"] = "WFI_IMAGE"
    rn_ref_model["meta"]["exposure"]["frame_time"] = 666

    if randomize:
        rn_ref_model["data"] = ((rng.random(shape) * rnoise).astype(np.float32),)
    else:
        rn_ref_model["data"] = np.ones(shape).astype(np.float32) * rnoise

    # Create temporary dark reference file
    # shape needs to be 3D but does not matter because the ramp fitting
    # step only uses the 2-D dark slope component
    dark_ref_model = DarkRefModel.create_fake_data(shape=(1, *shape))
    dark_ref_model["meta"]["instrument"]["detector"] = "WFI01"
    dark_ref_model["meta"]["instrument"]["name"] = "WFI"
    dark_ref_model["meta"]["reftype"] = "DARK"
    dark_ref_model["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    dark_ref_model["meta"]["exposure"]["type"] = "WFI_IMAGE"
    dark_ref_model["meta"]["exposure"]["frame_time"] = 666

    dark_ref_model["dark_slope"] = np.zeros(shape).astype(np.float32) * 0.01

    # return gainfile, readnoisefile
    return gain_ref_model, rn_ref_model, dark_ref_model
