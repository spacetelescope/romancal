from itertools import cycle

import pytest
import numpy as np

from jwst.datamodels import RampModel

from romancal.datamodels import GainModel, ReadNoiseModel
from romancal.jump import JumpStep

MAXIMUM_CORES = ['none', 'quarter', 'half', 'all']


@pytest.fixture(scope="module")
def generate_wfi_reffiles(tmpdir_factory):
    gainfile = str(tmpdir_factory.mktemp("ndata").join("gain.fits"))
    readnoisefile = str(tmpdir_factory.mktemp("ndata").join('readnoise.fits'))

    ingain = 6
    xsize = 20
    ysize = 20
    gain = np.ones(shape=(ysize, xsize), dtype=np.float64) * ingain
    gain_model = GainModel(data=gain)
    gain_model.meta.instrument.name = "WFI"

    inreadnoise = 5
    rnoise = np.ones(shape=(ysize, xsize), dtype=np.float64) * inreadnoise

    readnoise_model = ReadNoiseModel(data=rnoise)
    readnoise_model.meta.instrument.name = "WFI"
    readnoise_model.save(readnoisefile)

    return gainfile, readnoisefile


@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nints=1, nrows=1024, ncols=1032,
               nframes=1, grouptime=1.0, gain=1, deltatime=1):
        times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
        gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
        err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint32)

        rampmodel = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
        rampmodel.meta.instrument.name = 'WFI'
        rampmodel.meta.instrument.detector = 'WFI_MAGE'
        rampmodel.meta.instrument.filter = 'F480M'
        rampmodel.meta.observation.date = '2015-10-13'
        rampmodel.meta.exposure.type = 'WFI_IMAGE'
        rampmodel.meta.exposure.group_time = deltatime

        rampmodel.meta.exposure.frame_time = deltatime
        rampmodel.meta.exposure.ngroups = ngroups
        rampmodel.meta.exposure.group_time = deltatime
        rampmodel.meta.exposure.nframes = 1
        rampmodel.meta.exposure.groupgap = 0

        gain = GainModel(data=gain)
        gain.meta.instrument.name = 'WFI'

        rnmodel = ReadNoiseModel(data=read_noise)
        rnmodel.meta.instrument.name = 'WFI'

        return rampmodel, gdq, rnmodel, pixdq, err, gain

    return _setup


@pytest.mark.skip("Will be enabled when step can be called")
@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_wfi(generate_wfi_reffiles, setup_inputs, max_cores):
    override_gain, override_readnoise = generate_wfi_reffiles

    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = 7
    ngroups = 100
    CR_fraction = 5
    nrows = 20
    ncols = 20
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=nrows, ncols=ncols,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1, 89) if x % 5 == 0]
    CR_locs = [x for x in range(nrows * ncols) if x % CR_fraction == 0]
    CR_x_locs = [x % ncols for x in CR_locs]
    CR_y_locs = [int(x / nrows) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    # print("number of CRs {}".format(len(CR_x_locs)))

    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores=max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))
