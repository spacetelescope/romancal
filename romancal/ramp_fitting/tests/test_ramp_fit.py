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

RNG = np.random.default_rng(619)

MAXIMUM_CORES = ["none", "quarter", "half", "all"]

DO_NOT_USE = dqflags.group["DO_NOT_USE"]
JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]

dqflags = {
    "DO_NOT_USE": 1,
    "SATURATED": 2,
    "JUMP_DET": 4,
}


def generate_ramp_model(shape, deltatime=1):
    data = u.Quantity(
        (RNG.uniform(size=shape) * 0.5).astype(np.float32), u.DN, dtype=np.float32
    )
    err = u.Quantity(
        (RNG.uniform(size=shape) * 0.0001).astype(np.float32), u.DN, dtype=np.float32
    )
    pixdq = np.zeros(shape=shape[1:], dtype=np.uint32)
    gdq = np.zeros(shape=shape, dtype=np.uint8)

    dm_ramp = maker_utils.mk_ramp(shape=shape)
    dm_ramp.data = u.Quantity(data, u.DN, dtype=np.float32)
    dm_ramp.pixeldq = pixdq
    dm_ramp.groupdq = gdq
    dm_ramp.err = u.Quantity(err, u.DN, dtype=np.float32)

    dm_ramp.meta.exposure.frame_time = deltatime
    dm_ramp.meta.exposure.group_time = deltatime
    dm_ramp.meta.exposure.ngroups = shape[0]
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    ramp_model = RampModel(dm_ramp)

    return ramp_model


def generate_wfi_reffiles(shape, ingain=6):
    # Create temporary gain reference file
    gain_ref = maker_utils.mk_gain(shape=shape)

    gain_ref["meta"]["instrument"]["detector"] = "WFI01"
    gain_ref["meta"]["instrument"]["name"] = "WFI"
    gain_ref["meta"]["reftype"] = "GAIN"
    gain_ref["meta"]["useafter"] = Time("2022-01-01T11:11:11.111")

    gain_ref["data"] = u.Quantity(
        (RNG.uniform(size=shape) * 0.5).astype(np.float32) * ingain,
        u.electron / u.DN,
        dtype=np.float32,
    )
    gain_ref["dq"] = np.zeros(shape, dtype=np.uint16)
    gain_ref["err"] = u.Quantity(
        (RNG.uniform(size=shape) * 0.05).astype(np.float32),
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

    rn_ref["data"] = u.Quantity(
        (RNG.uniform(size=shape) * 0.01).astype(np.float32), u.DN, dtype=np.float32
    )

    rn_ref_model = ReadnoiseRefModel(rn_ref)

    # return gainfile, readnoisefile
    return gain_ref_model, rn_ref_model


@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason=(
        "Roman CRDS servers are not currently available outside the internal network"
    ),
)
@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_one_group_small_buffer_fit_ols(max_cores):
    ingain = 1.0
    deltatime = 1
    ngroups = 1
    xsize = 20
    ysize = 20
    shape = (ngroups, xsize, ysize)

    override_gain, override_readnoise = generate_wfi_reffiles(shape[1:], ingain)

    model1 = generate_ramp_model(shape, deltatime)

    model1.data[0, 15, 10] = 10.0 * model1.data.unit  # add single CR

    out_model = RampFitStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores=max_cores,
    )

    data = out_model.data.value

    # Index changes due to trimming of reference pixels
    np.testing.assert_allclose(data[11, 6], 10, 1e-6)


@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason=(
        "Roman CRDS servers are not currently available outside the internal network"
    ),
)
def test_multicore_ramp_fit_match():
    ingain = 1.0
    deltatime = 1
    ngroups = 4
    xsize = 20
    ysize = 20
    shape = (ngroups, xsize, ysize)

    override_gain, override_readnoise = generate_wfi_reffiles(shape[1:], ingain)

    model1 = generate_ramp_model(shape, deltatime)

    # gain or read noise are also modified in place in an important way (!)
    # so we make copies here so that we can get agreement.
    out_model = RampFitStep.call(
        model1.copy(),  # model1 is modified in place now.
        override_gain=override_gain.copy(),
        override_readnoise=override_readnoise.copy(),
        maximum_cores="none",
    )

    all_out_model = RampFitStep.call(
        model1.copy(),  # model1 is modified in place now.
        override_gain=override_gain.copy(),
        override_readnoise=override_readnoise.copy(),
        maximum_cores="all",
    )

    # Original ramp parameters
    np.testing.assert_allclose(out_model.data, all_out_model.data, 1e-6)
    np.testing.assert_allclose(out_model.err, all_out_model.err, 1e-6)
    np.testing.assert_allclose(out_model.amp33, all_out_model.amp33, 1e-6)
    np.testing.assert_allclose(
        out_model.border_ref_pix_left, all_out_model.border_ref_pix_left, 1e-6
    )
    np.testing.assert_allclose(
        out_model.border_ref_pix_right, all_out_model.border_ref_pix_right, 1e-6
    )
    np.testing.assert_allclose(
        out_model.border_ref_pix_top, all_out_model.border_ref_pix_top, 1e-6
    )
    np.testing.assert_allclose(
        out_model.border_ref_pix_bottom, all_out_model.border_ref_pix_bottom, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_left, all_out_model.dq_border_ref_pix_left, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_right, all_out_model.dq_border_ref_pix_right, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_top, all_out_model.dq_border_ref_pix_top, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_bottom, all_out_model.dq_border_ref_pix_bottom, 1e-6
    )

    # New rampfit parameters
    np.testing.assert_allclose(out_model.var_poisson, all_out_model.var_poisson, 1e-6)
    np.testing.assert_allclose(out_model.var_rnoise, all_out_model.var_rnoise, 1e-6)


@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason=(
        "Roman CRDS servers are not currently available outside the internal network"
    ),
)
@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_saturated_ramp_fit(max_cores):
    ingain = 1.0
    deltatime = 1
    ngroups = 4
    xsize = 20
    ysize = 20
    shape = (ngroups, xsize, ysize)

    # Create input model
    override_gain, override_readnoise = generate_wfi_reffiles(shape[1:], ingain)
    model1 = generate_ramp_model(shape, deltatime)

    # Set saturated flag
    model1.groupdq = model1.groupdq | SATURATED

    # Run ramp fit step
    out_model = RampFitStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores=max_cores,
    )

    # Test data and error arrays are zeroed out
    np.testing.assert_array_equal(out_model.data.value, 0)
    np.testing.assert_array_equal(out_model.err.value, 0)
    np.testing.assert_array_equal(out_model.var_poisson.value, 0)
    np.testing.assert_array_equal(out_model.var_rnoise.value, 0)

    # Test that all pixels are flagged saturated
    assert np.all(np.bitwise_and(out_model.dq, SATURATED) == SATURATED)

    # Test that original ramp parameters preserved
    np.testing.assert_allclose(out_model.amp33, model1.amp33, 1e-6)
    np.testing.assert_allclose(
        out_model.border_ref_pix_left, model1.border_ref_pix_left, 1e-6
    )
    np.testing.assert_allclose(
        out_model.border_ref_pix_right, model1.border_ref_pix_right, 1e-6
    )
    np.testing.assert_allclose(
        out_model.border_ref_pix_top, model1.border_ref_pix_top, 1e-6
    )
    np.testing.assert_allclose(
        out_model.border_ref_pix_bottom, model1.border_ref_pix_bottom, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_left, model1.dq_border_ref_pix_left, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_right, model1.dq_border_ref_pix_right, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_top, model1.dq_border_ref_pix_top, 1e-6
    )
    np.testing.assert_allclose(
        out_model.dq_border_ref_pix_bottom, model1.dq_border_ref_pix_bottom, 1e-6
    )

    # Test that an Image model was returned.
    assert type(out_model) == ImageModel

    # Test that the ramp fit step was labeled complete
    assert out_model.meta.cal_step.ramp_fit == "COMPLETE"
