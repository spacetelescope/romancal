"""Ramp Fitting tests involving MultiAccum Tables"""
import os

import numpy as np
import pytest
from astropy import units as u
from astropy.time import Time
from roman_datamodels import maker_utils
from roman_datamodels.datamodels import (
    GainRefModel,
    ImageModel,
    RampModel,
    ReadnoiseRefModel,
)

from romancal.lib import dqflags
from romancal.ramp_fitting import RampFitStep

DO_NOT_USE = dqflags.group["DO_NOT_USE"]
JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]

dqflags = {
    "DO_NOT_USE": 1,
    "SATURATED": 2,
    "JUMP_DET": 4,
}

# Basic resultant
#
# The read pattern is `[[1], [2], [3], [4]]`
# The total expected counts is 7.
# The resultants were generated with `romanisim.l1.apportion_counts_to_resultants(counts, read_pattern)`.
SIMPLE_RESULTANTS = np.array(
    [[[2., 2.],
      [5., 1.]],
     [[4., 5.],
      [6., 2.]],
     [[5., 6.],
      [7., 6.]],
     [[7., 7.],
      [7., 7.]]], dtype=np.float32)
SIMPLE_SLOPES = np.array([[0.52631587, 0.52631587], [0.23026317, 0.7236843 ]], dtype=np.float32)
GAIN_SLOPES = np.array([[2.667965, 2.667965], [1.133122, 3.50926 ]], dtype=np.float32)


@pytest.mark.parametrize(
    'resultants, ingain, rnoise, randomize, expected',
    [
        pytest.param(SIMPLE_RESULTANTS, 1, 0.01, False, SIMPLE_SLOPES, id='default'),     # No gain or noise
        pytest.param(SIMPLE_RESULTANTS, 5, 0.01, False, GAIN_SLOPES, id='extragain'),     # Increase gain
        pytest.param(SIMPLE_RESULTANTS, 1, 100., False, SIMPLE_SLOPES, id='extranoise'),  # Increase noise
    ],
)
def test_fit(resultants, ingain, rnoise, randomize, expected):
    """Test ramp fits"""
    ramp_model, gain_model, readnoise_model = make_data(resultants, ingain, rnoise, randomize)
    out_model = RampFitStep.call(
        ramp_model,
        algorithm='ols_cas21',
        override_gain=gain_model,
        override_readnoise=readnoise_model,
    )

    # Test for expectation
    data = out_model.data.value
    np.testing.assert_allclose(data, expected, 1e-6)


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
    gain_model, readnoise_model = generate_wfi_reffiles(ramp_model.shape[1:], ingain=ingain, rnoise=rnoise, randomize=randomize)

    return ramp_model, gain_model, readnoise_model


def model_from_resultants(resultants, read_pattern=None):
    """Create a RampModel from resultants

    Parameters
    ----------
    resultants : numpy.array.shape(reads, xdim, ydim)
        The resultants to fit.

    read_pattern : [[int[,...]][,...]]
        The read patter used to produce the resultants.
        If None, presume a basic read pattern
    """
    if read_pattern is None:
        read_pattern = [[idx + 1] for idx in range(resultants.shape[0])]

    # Full WFI image has reference pixels all around. Add those on.
    nrefpixs = 4
    full_wfi = np.ones((resultants.shape[0],
                        resultants.shape[1] + (nrefpixs * 2),
                        resultants.shape[2] + (nrefpixs * 2)),
                       dtype=np.float32)
    full_wfi[:,nrefpixs:-nrefpixs, nrefpixs:-nrefpixs] = resultants
    shape = full_wfi.shape

    pixdq = np.zeros(shape=shape[1:], dtype=np.uint32)
    err = np.zeros(shape=shape, dtype=np.float32)
    gdq = np.zeros(shape=shape, dtype=np.uint8)

    dm_ramp = maker_utils.mk_ramp(shape=shape)
    dm_ramp.data = u.Quantity(full_wfi, u.DN, dtype=np.float32)
    dm_ramp.pixeldq = pixdq
    dm_ramp.groupdq = gdq
    dm_ramp.err = u.Quantity(err, u.DN, dtype=np.float32)

    dm_ramp.meta.exposure.frame_time = 1
    dm_ramp.meta.exposure.ngroups = shape[0]
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    dm_ramp.meta.exposure.read_pattern = read_pattern

    ramp_model = RampModel(dm_ramp)

    return ramp_model


def generate_wfi_reffiles(shape, ingain=6, rnoise=0.01, randomize=True):
    """Create GainRefModel and ReadnoiseRefModel

    Parameters
    ----------
    shape : tuple
        Shape of the arrays

    ingain : float
        Maximum gain.

    rnoise : flota
        Maximum noise

    randomize : bool
        Randomize the gain and read noise data.
    """
    # Create temporary gain reference file
    gain_ref = maker_utils.mk_gain(shape=shape)

    gain_ref["meta"]["instrument"]["detector"] = "WFI01"
    gain_ref["meta"]["instrument"]["name"] = "WFI"
    gain_ref["meta"]["reftype"] = "GAIN"
    gain_ref["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    if randomize:
        gain_ref["data"] = u.Quantity(
            (np.random.random(shape) * 0.5).astype(np.float32) * ingain,
            u.electron / u.DN,
            dtype=np.float32,
        )
    else:
        gain_ref["data"] = u.Quantity(
            np.ones(shape).astype(np.float32) * ingain,
            u.electron / u.DN,
            dtype=np.float32,
        )
    gain_ref["dq"] = np.zeros(shape, dtype=np.uint16)
    gain_ref["err"] = u.Quantity(
        (np.random.random(shape) * 0.05).astype(np.float32),
        u.electron / u.DN,
        dtype=np.float32,
    )

    gain_ref_model = GainRefModel(gain_ref)

    # Create temporary readnoise reference file
    rn_ref = maker_utils.mk_readnoise(shape=shape)
    rn_ref["meta"]["instrument"]["detector"] = "WFI01"
    rn_ref["meta"]["instrument"]["name"] = "WFI"
    rn_ref["meta"]["reftype"] = "READNOISE"
    rn_ref["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    rn_ref["meta"]["exposure"]["type"] = "WFI_IMAGE"
    rn_ref["meta"]["exposure"]["frame_time"] = 666

    if randomize:
        rn_ref["data"] = u.Quantity(
            (np.random.random(shape) * rnoise).astype(np.float32), u.DN, dtype=np.float32
        )
    else:
        rn_ref["data"] = u.Quantity(
            np.ones(shape).astype(np.float32) * rnoise, u.DN, dtype=np.float32
        )

    rn_ref_model = ReadnoiseRefModel(rn_ref)

    # return gainfile, readnoisefile
    return gain_ref_model, rn_ref_model
