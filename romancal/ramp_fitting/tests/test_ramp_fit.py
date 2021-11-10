import pytest
import numpy as np

from astropy.time import Time

import roman_datamodels.stnode as rds
from roman_datamodels.datamodels import RampModel
from roman_datamodels.datamodels import GainRefModel
from roman_datamodels.datamodels import ReadnoiseRefModel
from roman_datamodels.testing import utils as testutil

from romancal.ramp_fitting import RampFitStep
from romancal.lib import dqflags


# MAXIMUM_CORES = ['none', 'quarter', 'half', 'all']
MAXIMUM_CORES = ['none']  # initial testing only

DO_NOT_USE = dqflags.group['DO_NOT_USE']
JUMP_DET = dqflags.group['JUMP_DET']
SATURATED = dqflags.group['SATURATED']

dqflags = {
    "DO_NOT_USE": 1,
    "SATURATED": 2,
    "JUMP_DET": 4,
}


@pytest.fixture
def setup_inputs():

    def _setup(ngroups=10, nrows=20, ncols=20, deltatime=1):

        data = np.zeros(shape=(ngroups, nrows, ncols), dtype=np.float32)
        err = np.ones(shape=(nrows, ncols), dtype=np.float32)
        pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
        gdq = np.zeros(shape=(ngroups, nrows, ncols), dtype=np.uint8)
        csize = (ngroups, nrows, ncols)

        dm_ramp = testutil.mk_ramp(csize)
        dm_ramp.data = data
        dm_ramp.pixeldq = pixdq
        dm_ramp.groupdq = gdq
        dm_ramp.err = err

        dm_ramp.meta.exposure.frame_time = deltatime
        dm_ramp.meta.exposure.ngroups = ngroups
        dm_ramp.meta.exposure.nframes = 1
        dm_ramp.meta.exposure.groupgap = 0

        return dm_ramp

    return _setup


@pytest.fixture(scope="module")
def generate_wfi_reffiles(tmpdir_factory):

    gainfile = str(tmpdir_factory.mktemp("ndata").join("gain.asdf"))
    readnoisefile = str(tmpdir_factory.mktemp("ndata").join('readnoise.asdf'))

    ingain = 6
    xsize = 20
    ysize = 20

    shape = (ysize, xsize)

    # Create temporary gain reference file
    gain_ref = rds.GainRef()

    meta = {}
    testutil.add_ref_common(meta)
    meta['instrument']['detector'] = 'WFI01'
    meta['instrument']['name'] = 'WF1'
    meta['author'] = 'John Doe'
    meta['reftype'] = 'GAIN'
    meta['pedigree'] = 'DUMMY'
    meta['description'] = 'DUMMY'
    meta['useafter'] = Time('2022-01-01T11:11:11.111')

    gain_ref['meta'] = meta
    gain_ref['data'] = np.ones(shape, dtype=np.float32) * ingain
    gain_ref['dq'] = np.zeros(shape, dtype=np.uint16)
    gain_ref['err'] = (np.random.random(shape) * 0.05).astype(np.float64)

    gain_ref_model = GainRefModel(gain_ref)
    gain_ref_model.save(gainfile)

    # Create temporary readnoise reference file
    rn_ref = rds.ReadnoiseRef()
    meta = {}
    testutil.add_ref_common(meta)
    meta['instrument']['detector'] = 'WFI01'
    meta['instrument']['name'] = 'WF1'
    meta['author'] = 'John Doe'
    meta['reftype'] = 'READNOISE'
    meta['pedigree'] = 'DUMMY'
    meta['description'] = 'DUMMY'
    meta['useafter'] = Time('2022-01-01T11:11:11.111')

    exposure = {}
    exposure['type'] = 'WFI_IMAGE'
    exposure['frame_time'] = 666

    rn_ref['meta'] = meta
    rn_ref['meta']['exposure'] = exposure
    rn_ref['data'] = np.ones(shape, dtype=np.float32)

    rn_ref_model = ReadnoiseRefModel(rn_ref)
    rn_ref_model.save(readnoisefile)

    return gainfile, readnoisefile


@pytest.fixture(scope="module")
def generate_wfi_inputfile(tmpdir_factory):

    infile = str(tmpdir_factory.mktemp("ndata").join("input.asdf"))

    return infile


# @pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
@pytest.mark.xfail
def test_one_group_small_buffer_fit_ols(generate_wfi_reffiles,
                                        generate_wfi_inputfile, max_cores,
                                        setup_inputs):

    override_gain, override_readnoise = generate_wfi_reffiles
    inputfile = generate_wfi_inputfile

    grouptime = 1.
    ingain = 1.
    inreadnoise = 10.
    ngroups = 1
    xsize = 20
    ysize = 20

    model1 = setup_inputs(ngroups=ngroups, nrows=ysize, ncols=xsize,
                          gain=ingain, readnoise=inreadnoise,
                          deltatime=grouptime)

    model1.data[0, 15, 10] = 10.0  # add single CR

    a_ramp_model = RampModel(model1)
    a_ramp_model.save(inputfile)

    out_model, int_model = \
        RampFitStep.call(inputfile, override_gain=override_gain,
                         override_readnoise=override_readnoise,
                         maximum_cores=max_cores)

    data = out_model.data

    np.testing.assert_allclose(data[15, 10], 10.0, 1e-6)
