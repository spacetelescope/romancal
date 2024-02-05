"""
 Unit tests for the Roman jump step code
"""

from itertools import cycle

import numpy as np
import pytest
import roman_datamodels.stnode as rds
from astropy import units as u
from astropy.time import Time
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils
from roman_datamodels.datamodels import GainRefModel, ReadnoiseRefModel

from romancal.jump import JumpStep

RNG = np.random.default_rng(42)

MAXIMUM_CORES = ["none", "quarter", "half", "all"]


@pytest.fixture(scope="module")
def generate_wfi_reffiles(tmpdir_factory):
    gainfile = str(tmpdir_factory.mktemp("ndata").join("gain.asdf"))
    readnoisefile = str(tmpdir_factory.mktemp("ndata").join("readnoise.asdf"))

    ingain = 6
    xsize = 20
    ysize = 20

    shape = (ysize, xsize)

    # Create temporary gain reference file
    gain_ref = rds.GainRef()
    meta = maker_utils.mk_ref_common("GAIN")
    meta["instrument"]["detector"] = "WFI01"
    meta["instrument"]["name"] = "WFI"
    meta["author"] = "John Doe"
    meta["pedigree"] = "DUMMY"
    meta["description"] = "DUMMY"
    meta["useafter"] = Time("2022-01-01T11:11:11.111")

    gain_ref["meta"] = meta
    gain_ref["data"] = u.Quantity(
        np.ones(shape, dtype=np.float32) * ingain, u.electron / u.DN, dtype=np.float32
    )
    gain_ref["dq"] = np.zeros(shape, dtype=np.uint16)
    gain_ref["err"] = u.Quantity(
        (RNG.uniform(size=shape) * 0.05).astype(np.float64),
        u.electron / u.DN,
        dtype=np.float64,
    )

    gain_ref_model = GainRefModel(gain_ref)
    gain_ref_model.save(gainfile)

    # Create temporary readnoise reference file
    rn_ref = rds.ReadnoiseRef()
    meta = maker_utils.mk_ref_common("READNOISE")
    meta["instrument"]["detector"] = "WFI01"
    meta["instrument"]["name"] = "WFI"
    meta["author"] = "John Doe"
    meta["pedigree"] = "DUMMY"
    meta["description"] = "DUMMY"
    meta["useafter"] = Time("2022-01-01T11:11:11.111")
    meta["exposure"] = {}
    meta["exposure"]["type"] = "WFI_IMAGE"
    meta["exposure"]["p_exptype"] = "WFI_IMAGE|WFI_GRISM|WFI_PRISM|"

    rn_ref["meta"] = meta
    rn_ref["data"] = u.Quantity(
        np.ones(shape, dtype=np.float32), u.DN, dtype=np.float32
    )
    rn_ref["dq"] = np.zeros(shape, dtype=np.uint16)
    rn_ref["err"] = u.Quantity(
        (RNG.uniform(size=shape) * 0.05).astype(np.float64), u.DN, dtype=np.float64
    )

    rn_ref_model = ReadnoiseRefModel(rn_ref)
    rn_ref_model.save(readnoisefile)

    return gainfile, readnoisefile


@pytest.fixture
def setup_inputs():
    def _setup(
        ngroups=10,
        readnoise=10,
        nrows=20,
        ncols=20,
        nframes=1,
        grouptime=1.0,
        gain=1,
        deltatime=1,
    ):
        err = np.ones(shape=(ngroups, nrows, ncols), dtype=np.float32)
        data = np.zeros(shape=(ngroups, nrows, ncols), dtype=np.float32)
        gdq = np.zeros(shape=(ngroups, nrows, ncols), dtype=np.uint8)
        pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)

        csize = (ngroups, nrows, ncols)
        dm_ramp = rdm.RampModel(maker_utils.mk_ramp(shape=csize))

        dm_ramp.meta.instrument.name = "WFI"
        dm_ramp.meta.instrument.optical_element = "F158"

        dm_ramp.data = u.Quantity(data + 6.0, u.DN, dtype=np.float32)
        dm_ramp.pixeldq = pixdq
        dm_ramp.groupdq = gdq
        dm_ramp.err = u.Quantity(err, u.DN, dtype=np.float32)

        dm_ramp.meta.exposure.type = "WFI_IMAGE"
        dm_ramp.meta.exposure.group_time = deltatime
        dm_ramp.meta.exposure.frame_time = deltatime
        dm_ramp.meta.exposure.ngroups = ngroups
        dm_ramp.meta.exposure.nframes = 1
        dm_ramp.meta.exposure.groupgap = 0
        dm_ramp.meta.cal_step["dq_init"] = "INCOMPLETE"
        dm_ramp.meta.cal_step["jump"] = "INCOMPLETE"

        return dm_ramp

    return _setup


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_one_CR(generate_wfi_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_wfi_reffiles

    grouptime = 3.0
    deltaDN = 5
    ingain = 200
    inreadnoise = 7.0
    ngroups = 100
    CR_fraction = 3
    xsize = 20
    ysize = 20

    model1 = setup_inputs(
        ngroups=ngroups,
        nrows=ysize,
        ncols=xsize,
        gain=ingain,
        readnoise=inreadnoise,
        deltatime=grouptime,
    )

    for i in range(ngroups):
        model1.data[i, :, :] = deltaDN * i * model1.data.unit

    first_CR_group_locs = [x for x in range(1, 7) if x % 5 == 0]

    CR_locs = [x for x in range(xsize * ysize) if x % CR_fraction == 0]
    CR_x_locs = [x % ysize for x in CR_locs]
    CR_y_locs = [int(x / xsize) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)

    # Add CRs to the SCI data
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[CR_group:, CR_y_locs[i], CR_x_locs[i]] = (
            model1.data[CR_group:, CR_y_locs[i], CR_x_locs[i]]
            + 500.0 * model1.data.unit
        )

    out_model = JumpStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores=max_cores,
        use_ramp_jump_detection=False,
    )

    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)

        assert 4 == np.max(out_model.groupdq[CR_group, CR_y_locs[i], CR_x_locs[i]])


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_two_CRs(generate_wfi_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_wfi_reffiles

    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = 7.0
    ngroups = 100
    CR_fraction = 5
    xsize = 20
    ysize = 20

    model1 = setup_inputs(
        ngroups=ngroups,
        nrows=ysize,
        ncols=xsize,
        gain=ingain,
        readnoise=inreadnoise,
        deltatime=grouptime,
    )

    for i in range(ngroups):
        model1.data[i, :, :] = deltaDN * i * model1.data.unit

    first_CR_group_locs = [x for x in range(1, 7) if x % 5 == 0]
    CR_locs = [x for x in range(xsize * ysize) if x % CR_fraction == 0]
    CR_x_locs = [x % ysize for x in CR_locs]
    CR_y_locs = [int(x / xsize) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)

    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)

        model1.data[CR_group:, CR_y_locs[i], CR_x_locs[i]] = (
            model1.data[CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500 * model1.data.unit
        )
        model1.data[CR_group + 8 :, CR_y_locs[i], CR_x_locs[i]] = (
            model1.data[CR_group + 8 :, CR_y_locs[i], CR_x_locs[i]]
            + 700 * model1.data.unit
        )

    out_model = JumpStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores=max_cores,
        use_ramp_jump_detection=False,
    )

    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)

        assert 4 == np.max(out_model.groupdq[CR_group, CR_y_locs[i], CR_x_locs[i]])
        assert 4 == np.max(out_model.groupdq[CR_group + 8, CR_y_locs[i], CR_x_locs[i]])


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_two_group_integration(generate_wfi_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_wfi_reffiles
    grouptime = 3.0
    ingain = 6
    inreadnoise = 7
    ngroups = 2
    xsize = 20
    ysize = 20
    model1 = setup_inputs(
        ngroups=ngroups,
        nrows=ysize,
        ncols=xsize,
        gain=ingain,
        readnoise=inreadnoise,
        deltatime=grouptime,
    )

    out_model = JumpStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores=max_cores,
    )

    assert out_model.meta.cal_step["jump"] == "SKIPPED"


@pytest.mark.parametrize("use_jump", [False, True])
def test_four_group_integration(generate_wfi_reffiles, setup_inputs, use_jump):
    override_gain, override_readnoise = generate_wfi_reffiles
    grouptime = 3.0
    ingain = 6
    inreadnoise = 7
    ngroups = 4
    xsize = 20
    ysize = 20
    model1 = setup_inputs(
        ngroups=ngroups,
        nrows=ysize,
        ncols=xsize,
        gain=ingain,
        readnoise=inreadnoise,
        deltatime=grouptime,
    )

    out_model = JumpStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores="none",
        use_ramp_jump_detection=use_jump,
    )

    if use_jump:
        assert out_model.meta.cal_step["jump"] == "SKIPPED"
    else:
        assert out_model.meta.cal_step["jump"] == "COMPLETE"


@pytest.mark.parametrize("use_jump", [False, True])
def test_three_group_integration(generate_wfi_reffiles, setup_inputs, use_jump):
    override_gain, override_readnoise = generate_wfi_reffiles
    grouptime = 3.0
    ingain = 6
    inreadnoise = 7
    ngroups = 3
    xsize = 20
    ysize = 20
    model1 = setup_inputs(
        ngroups=ngroups,
        nrows=ysize,
        ncols=xsize,
        gain=ingain,
        readnoise=inreadnoise,
        deltatime=grouptime,
    )

    out_model = JumpStep.call(
        model1,
        override_gain=override_gain,
        override_readnoise=override_readnoise,
        maximum_cores="none",
        use_ramp_jump_detection=use_jump,
    )

    if use_jump:
        assert out_model.meta.cal_step.jump == "SKIPPED"
    else:
        assert out_model.meta.cal_step.jump == "COMPLETE"
