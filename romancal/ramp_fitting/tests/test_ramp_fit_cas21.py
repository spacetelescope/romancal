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


@pytest.mark.skipif(
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
        ingain, deltatime, ngroups, xsize, ysize = request.param
    else:
        ingain = 1
        deltatime = 1
        ngroups = 4
        xsize = 20
        ysize = 20
    shape = (ngroups, xsize, ysize)

    image = generate_ramp_model(shape, deltatime)
    gain, readnoise = generate_wfi_reffiles(shape[1:], ingain)

    return image, gain, readnoise


# #########
# Utilities
# #########

def generate_ramp_model(shape, deltatime=1):
    data = (np.random.random(shape) * 0.5).astype(np.float32)
    err = (np.random.random(shape) * 0.0001).astype(np.float32)
    pixdq = np.zeros(shape=shape[1:], dtype=np.uint32)
    gdq = np.zeros(shape=shape, dtype=np.uint8)

    dm_ramp = maker_utils.mk_ramp(shape)
    dm_ramp.data = u.Quantity(data, u.DN, dtype=np.float32)
    dm_ramp.data = u.Quantity(data, u.DN, dtype=np.float32)
    dm_ramp.pixeldq = pixdq
    dm_ramp.groupdq = gdq
    dm_ramp.err = u.Quantity(err, u.DN, dtype=np.float32)

    dm_ramp.meta.exposure.frame_time = deltatime
    dm_ramp.meta.exposure.ngroups = shape[0]
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    ramp_model = RampModel(dm_ramp)

    return ramp_model


def generate_wfi_reffiles(shape, ingain=6):
    # Create temporary gain reference file
    gain_ref = maker_utils.mk_gain(shape)

    gain_ref["meta"]["instrument"]["detector"] = "WFI01"
    gain_ref["meta"]["instrument"]["name"] = "WFI"
    gain_ref["meta"]["reftype"] = "GAIN"
    gain_ref["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    gain_ref["data"] = u.Quantity(
        (np.random.random(shape) * 0.5).astype(np.float32) * ingain,
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
    rn_ref = maker_utils.mk_readnoise(shape)
    rn_ref["meta"]["instrument"]["detector"] = "WFI01"
    rn_ref["meta"]["instrument"]["name"] = "WFI"
    rn_ref["meta"]["reftype"] = "READNOISE"
    rn_ref["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    rn_ref["meta"]["exposure"]["type"] = "WFI_IMAGE"
    rn_ref["meta"]["exposure"]["frame_time"] = 666

    rn_ref["data"] = u.Quantity(
        (np.random.random(shape) * 0.01).astype(np.float32), u.DN, dtype=np.float32
    )

    rn_ref_model = ReadnoiseRefModel(rn_ref)

    # return gainfile, readnoisefile
    return gain_ref_model, rn_ref_model
