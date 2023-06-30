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

MAXIMUM_CORES = ["none", "quarter", "half", "all"]

DO_NOT_USE = dqflags.group["DO_NOT_USE"]
JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]

dqflags = {
    "DO_NOT_USE": 1,
    "SATURATED": 2,
    "JUMP_DET": 4,
}


pytestmark = pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason=(
        "Roman CRDS servers are not currently available outside the internal network"
    ),
)


def test_ols_cas21_default(make_data):
    model, override_gain, override_readnoise = make_data

    out_model = RampFitStep.call(
        model,
        algorithm='ols_cas21',
        override_gain=override_gain,
        override_readnoise=override_readnoise,
    )

    data = out_model.data.value

    # Index changes due to trimming of reference pixels
    np.testing.assert_allclose(data[11, 6], 10.0, 1e-6)


def test_ols_cas21_fixed(make_data_fixed):
    model, override_gain, override_readnoise, slopes = make_data_fixed

    out_model = RampFitStep.call(
        model,
        algorithm='ols_cas21',
        override_gain=override_gain,
        override_readnoise=override_readnoise,
    )

    # Test for expectation
    data = out_model.data.value
    np.testing.assert_allclose(data, slopes, 1e-6)


# ########
# fixtures
# ########
@pytest.fixture(scope='module')
def make_data(request):
    """Create test input data

    Parameters
    ----------
    request.param : (ingain, deltatime, ngroups, xsize, ysize)
        If specified, set the parameters of the created data.
        If not specified, defaults are used.

    Returns
    -------
    image, gain, readnoise : ImageModel, GainRefModel, ReadnoiseRefModel
        Input image and related references
    """
    if getattr(request, 'param', None):
        ingain, deltatime, ngroups, xsize, ysize, randomize = request.param
    else:
        ingain = 1
        deltatime = 1
        ngroups = 4
        xsize = 20
        ysize = 20
        randomize = True
    shape = (ngroups, xsize, ysize)

    image = generate_ramp_model(shape, deltatime)
    gain, readnoise = generate_wfi_reffiles(shape[1:], ingain=ingain, randomize=randomize)

    return image, gain, readnoise


@pytest.fixture(scope='module')
def make_data_fixed():
    """Create test input data

    Parameters
    ----------
    request.param : (ingain, deltatime, ngroups, xsize, ysize)
        If specified, set the parameters of the created data.
        If not specified, defaults are used.

    Returns
    -------
    image, gain, readnoise, expected : ImageModel, GainRefModel, ReadnoiseRefModel, numpy.array
        Input image, related references, and expected slopes
    """
    expected = np.array([[0.52631587, 0.52631587], [0.23026317, 0.7236843 ]], dtype=np.float32)

    image = generate_ramp_model_constant()
    gain, readnoise = generate_wfi_reffiles(image.shape[1:], ingain=1, randomize=False)

    return image, gain, readnoise, expected


# #########
# Utilities
# #########

def generate_ramp_model(shape, deltatime=1, read_pattern=None):
    data = (np.random.random(shape) * 0.5).astype(np.float32)
    err = (np.random.random(shape) * 0.0001).astype(np.float32)
    pixdq = np.zeros(shape=shape[1:], dtype=np.uint32)
    gdq = np.zeros(shape=shape, dtype=np.uint8)

    dm_ramp = maker_utils.mk_ramp(shape=shape)
    dm_ramp.data = u.Quantity(data, u.DN, dtype=np.float32)
    dm_ramp.pixeldq = pixdq
    dm_ramp.groupdq = gdq
    dm_ramp.err = u.Quantity(err, u.DN, dtype=np.float32)

    dm_ramp.meta.exposure.frame_time = deltatime
    dm_ramp.meta.exposure.ngroups = shape[0]
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # If not read pattern, presume just an evenly defined set.
    if read_pattern is None:
        read_pattern = [[idx] for idx in range(1, shape[0] + 1)]
    dm_ramp.meta.exposure.read_pattern = read_pattern

    ramp_model = RampModel(dm_ramp)

    return ramp_model


def generate_ramp_model_constant():
    """Generate simple ramp model

    The read pattern is `[[1], [2], [3], [4]]`
    The total expected counts is 7.
    The resultants were generated with `romanisim.l1.apportion_counts_to_resultants(counts, read_pattern)`.
    """
    resultants = np.array(
        [[[2., 2.],
          [5., 1.]],
         [[4., 5.],
          [6., 2.]],
         [[5., 6.],
          [7., 6.]],
         [[7., 7.],
          [7., 7.]]], dtype=np.float32)

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

    dm_ramp.meta.exposure.read_pattern = [[1], [2], [3], [4]]

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
