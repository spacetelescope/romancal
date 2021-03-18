"""
Unit tests for dark current correction
"""

import warnings
import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.time import Time

from romancal.dark_current.dark_sub import (
    average_dark_frames,
    do_correction as darkcorr
    )
from romancal.datamodels import RampModel, DarkModel, dqflags


# Define frame_time and number of groups in the generated dark reffile
TFRAME = 10.0
NGROUPS_DARK = 10


def _params():
    """ Returns list of tuples, one for several combinations of ngroups,
        nframes, and nskip (which are some of the parameter which will later
        be retrieved from the MA table), generating parameters for
        test_frame_averaging.  Parameters are the following:

        (ma_tab, ngroups, nframes, groupgap, nrows, ncols)
    """
    # Dictionary of pseudo readout patterns
    ma_tab_infos = dict(
        RP1=dict(ngroups=20, nframes=8, nskip=0),
        RP2=dict(ngroups=32, nframes=2, nskip=0),
        RP3=dict(ngroups=12, nframes=6, nskip=0),
    )

    params = []
    ngroups = 3

    # WFI is 4096x4096, but we reduce the size to 20x20 for speed/memory
    nrows = 20
    ncols = 20
    for ma_tab_info, values in ma_tab_infos.items():
        params.append((ma_tab_info, ngroups, values['nframes'],
                       values['nskip'], nrows, ncols))

    return params


@pytest.mark.parametrize('ma_tab_info, ngroups, nframes, groupgap, nrows,'
                         'ncols', _params())
def test_frame_averaging(setup_nrc_cube, ma_tab_info, ngroups, nframes,
                         groupgap, nrows, ncols):

    '''Check that if nframes>1 or groupgap>0, then the pipeline reconstructs
       the dark reference file to match the frame averaging and groupgap
       settings of the exposure.'''

    # Create data and dark model
    data, dark = setup_nrc_cube(ma_tab_info, ngroups, nframes, groupgap,
                                nrows, ncols)

    # Add ramp values to dark model data array
    dark.data[:, 10, 10] = np.arange(0, NGROUPS_DARK)
    dark.err[:, 10, 10] = np.arange(10, NGROUPS_DARK + 10)

    # Run the pipeline's averaging function
    avg_dark = average_dark_frames(dark, ngroups, nframes, groupgap)

    # Group input groups into collections of frames which will be averaged
    total_frames = (nframes * ngroups) + (groupgap * (ngroups-1))

    # Get starting/ending indexes of the input groups to be averaged
    gstrt_ind = np.arange(0, total_frames, nframes + groupgap)
    gend_ind = gstrt_ind + nframes

    # Prepare arrays to hold results of averaging
    manual_avg = np.zeros((ngroups), dtype=np.float32)
    manual_errs = np.zeros((ngroups), dtype=np.float32)

    # Manually average the input data to compare with pipeline output
    for newgp, gstart, gend in zip(range(ngroups), gstrt_ind, gend_ind):

        # Average the data frames
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            newframe = np.mean(dark.data[gstart:gend, 10, 10])

            manual_avg[newgp] = newframe

        # ERR arrays will be quadratic sum of error values
        err_array = np.sum(dark.err[gstart:gend, 10, 10]**2)
        manual_errs[newgp] = np.sqrt(err_array) / (gend - gstart)

    # Check that pipeline output matches manual averaging results
    assert_allclose(manual_avg, avg_dark.data[:, 10, 10], rtol=1e-5)
    assert_allclose(manual_errs, avg_dark.err[:, 10, 10], rtol=1e-5)

    # Check that meta data was properly updated
    assert avg_dark.meta.exposure.nframes == nframes
    assert avg_dark.meta.exposure.ngroups == ngroups
    assert avg_dark.meta.exposure.groupgap == groupgap


def test_more_sci_frames(make_rampmodel, make_darkmodel):
    '''Check that data is unchanged if there are more frames in the science
    data is than in the dark reference file and verify that when the dark is
    not applied, the data is correctly flagged as such'''

    # size of integration
    ngroups = 7
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(ngroups, ysize, xsize)

    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[0, i] = i

    refgroups = 5
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # test that the science data are not changed; input file = output file
    np.testing.assert_array_equal(dm_ramp.data, outfile.data)

    darkstatus = outfile.meta.cal_step.dark_sub
    assert darkstatus == 'SKIPPED'


def test_sub_by_frame(make_rampmodel, make_darkmodel):
    '''Check that if NFRAMES=1 and GROUPGAP=0 for the science data, the
    dark reference data are directly subtracted frame by frame'''

    # size of integration
    ngroups = 6
    xsize = 3
    ysize = 3

    # create raw input data for step
    dm_ramp = make_rampmodel(ngroups, ysize, xsize)

    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[i] = i

    refgroups = 8
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[i] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)
    outdata = outfile.data

    diff = dm_ramp.data - dark.data[:ngroups]

    # test that the output data file is equal to the difference found when
    #     subtracting ref file from sci file
    np.testing.assert_array_equal(outdata, diff,
                                  err_msg='dark file should be subtracted'
                                          ' from sci file ')


def test_nan(make_rampmodel, make_darkmodel):
    '''Verify that when a dark has NaNs, these are correctly
    assumed as zero and the PIXELDQ is set properly'''

    # size of integration
    ngroups = 10
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(ngroups, ysize, xsize)

    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[i, :, :] = i + 0.

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # set NaN in dark file
    dark.data[5, 100, 100] = np.nan

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # test that the NaN dark reference pixel was set to 0 (nothing subtracted)

    assert outfile.data[5, 100, 100] == 5.0


def test_dq_combine(make_rampmodel, make_darkmodel):
    '''Verify that the DQ array of the dark is correctly combined
    with the PIXELDQ array of the science data.'''

    # size of integration
    ngroups = 5
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(ngroups, ysize, xsize)

    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(1, ngroups-1):
        dm_ramp.data[i, :, :] = i

    # create dark reference file model with more frames than science data
    refgroups = 7
    dark = make_darkmodel(refgroups, ysize, xsize)

    jump_det = dqflags.pixel['JUMP_DET']
    saturated = dqflags.pixel['SATURATED']
    do_not_use = dqflags.pixel['DO_NOT_USE']

    # populate dq flags of sci pixeldq and reference dq
    dm_ramp.pixeldq[50, 50] = jump_det
    dm_ramp.pixeldq[50, 51] = saturated

    dark.dq[50, 50] = do_not_use
    dark.dq[50, 51] = do_not_use

    # run correction step
    outfile = darkcorr(dm_ramp, dark)

    # check that dq flags were correctly added
    assert outfile.pixeldq[50, 50] == np.bitwise_or(jump_det, do_not_use)
    assert outfile.pixeldq[50, 51] == np.bitwise_or(saturated, do_not_use)


def test_frame_avg(make_rampmodel, make_darkmodel):
    '''Check that if NFRAMES>1 or GROUPGAP>0, the frame-averaged dark data are
    subtracted group-by-group from science data groups and the ERR arrays are
    not modified'''

    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(ngroups, ysize, xsize)

    dm_ramp.meta.exposure.nframes = 4
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i] = i + 1

    # create dark reference file model
    refgroups = 20  # This needs to be 20 groups for the calculations to work
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # dark frames should be averaged in groups of 4 frames
    # this will result in average values of 0.75, 0.775, 0.8, and 0.825
    # these values are then subtracted from frame values of 1, 2, 3 and 4
    assert outfile.data[0, 0, 0] == pytest.approx(0.25)
    assert outfile.data[0, 1, 0] == pytest.approx(1.225)
    assert outfile.data[0, 2, 0] == pytest.approx(2.2)
    assert outfile.data[0, 3, 0] == pytest.approx(3.175)

    # check that the error array is not modified.
    np.testing.assert_array_equal(outfile.err[:, :], 0)


@pytest.fixture(scope='function')
def make_rampmodel():
    '''Make Ramp model for testing'''

    def _ramp(ngroups, ysize, xsize):

        # create the data and groupdq arrays
        csize = (ngroups, ysize, xsize)

        data = np.full(csize, 1.0)

        # create a Roman datamodel for WFI data
        dm_ramp = RampModel(data=data)
        dm_ramp.meta.instrument.name = 'WFI'
        dm_ramp.meta.observation.datee = Time('2015-10-01T00:00:00')
        dm_ramp.meta.observation.time = Time('2015-10-01T01:01:00')
        dm_ramp.meta.description = 'Fake data.'

        return dm_ramp

    return _ramp


@pytest.fixture(scope='function')
def make_darkmodel():
    '''Make WFI dark model for testing'''

    def _dark(ngroups, ysize, xsize):
        # create the data and groupdq arrays
        csize = (ngroups, ysize, xsize)
        data = np.full(csize, 1.0)

        # create a dark datamodel for WFI data
        dark = DarkModel(data=data)
        dark.meta.instrument.name = 'WFI'

        dark.meta.exposure.nframes = 1
        dark.meta.exposure.groupgap = 0
        dark.meta.description = 'Fake data.'
        dark.meta.reftype = 'DARK'
        dark.meta.author = 'Alicia'
        dark.meta.pedigree = 'Dummy'
        dark.meta.useafter = Time('2015-10-01T00:00:00')

        return dark

    return _dark


@pytest.fixture(scope='function')
def setup_nrc_cube():
    '''Set up fake NIRCam data to test.'''

    def _cube(ma_tab_info, ngroups, nframes, groupgap, nrows, ncols):
        data_model = RampModel((ngroups, nrows, ncols))

        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.exposure.groupgap = groupgap
        data_model.meta.exposure.nframes = nframes
        data_model.meta.exposure.frame_time = TFRAME
        data_model.meta.exposure.group_time = (nframes + groupgap) * TFRAME
        data_model.meta.instrument.name = 'WFI'
        data_model.meta.instrument.detector = 'WFI01'

        data_model.meta.observation.date = Time('2017-11-01T00:00:00')
        data_model.meta.observation.time = Time('2017-11-01T05:13:00')

        dark_model = DarkModel((NGROUPS_DARK, 2048, 2048))

        dark_model.meta.exposure.ngroups = NGROUPS_DARK
        dark_model.meta.exposure.groupgap = 0
        dark_model.meta.exposure.nframes = 1
        dark_model.meta.instrument.name = 'WFI'
        dark_model.meta.description = 'Fake data.'
        dark_model.meta.telescope = 'ROMAN'
        dark_model.meta.reftype = 'DARK'
        dark_model.meta.author = 'Alicia'
        dark_model.meta.pedigree = 'Dummy'
        dark_model.meta.useafter = Time('2016-10-01T00:00:00')

        return data_model, dark_model

    return _cube
