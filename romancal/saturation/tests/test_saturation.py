"""

Unit tests for saturation flagging

"""

import pytest
import numpy as np

from romancal.saturation.saturation import flag_saturation
from romancal.lib import dqflags
from roman_datamodels.testing import utils as testutil


def test_basic_saturation_flagging(setup_wfi_datamodels):
    '''Check that the saturation flag is set when a pixel value is above the
       threshold given by the reference file.'''

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000
    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 5, 5] = 0
    data.data[1, 5, 5] = 20000
    data.data[2, 5, 5] = 40000
    data.data[3, 5, 5] = 60000   # Signal reaches saturation limit
    data.data[4, 5, 5] = 62000

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Make sure that groups with signal > saturation limit get flagged
    satindex = np.argmax(output.data[:, 5, 5] == satvalue)
    assert np.all(output.groupdq[satindex:, 5, 5] == dqflags.group['SATURATED'])

def test_ad_floor_flagging(setup_wfi_datamodels):
    """Check that the ad_floor flag is set when a pixel value is zero or
    negative."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[1, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[2, 5, 5] = 20
    data.data[3, 5, 5] = 40
    data.data[4, 5, 5] = 60

    # frames that should be flagged as saturation (low)
    satindxs = [0, 1]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Check if the right frames are flagged as saturated
    assert np.all(output.groupdq[satindxs, 5, 5]
                  == dqflags.group['DO_NOT_USE'] | dqflags.group['AD_FLOOR'])


def test_ad_floor_and_saturation_flagging(setup_wfi_datamodels):
    """Check that the ad_floor flag is set when a pixel value is zero or
    negative and the saturation flag when the pixel is above the saturation threshold."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[1, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[2, 5, 5] = 20
    data.data[3, 5, 5] = 40
    data.data[4, 5, 5] = 61000  # Signal above the saturation threshold

    # frames that should be flagged as ad_floor
    floorindxs = [0, 1]
    # frames that should be flagged as saturation
    satindxs = [4]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Check if the right frames are flagged as ad_floor
    assert np.all(output.groupdq[floorindxs, 5, 5]
                  == dqflags.group['DO_NOT_USE'] | dqflags.group['AD_FLOOR'])
    # Check if the right frames are flagged as saturated
    assert np.all(output.groupdq[satindxs, 5, 5] == dqflags.group['SATURATED'])


def test_signal_fluctuation_flagging(setup_wfi_datamodels):
    '''Check that once a pixel is flagged as saturated in a group, all
       subsequent groups should also be flagged as saturated, even if the
       signal value drops back below saturation.'''

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 5, 5] = 10
    data.data[1, 5, 5] = 20000
    data.data[2, 5, 5] = 40000
    data.data[3, 5, 5] = 60000   # Signal reaches saturation limit
    data.data[4, 5, 5] = 40000   # Signal drops below saturation limit

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Make sure that all groups after first saturated group are flagged
    satindex = np.argmax(output.data[:, 5, 5] == satvalue)
    assert np.all(output.groupdq[satindex:, 5, 5] == dqflags.group['SATURATED'])


def test_all_groups_saturated(setup_wfi_datamodels):
    '''Check case where all groups are saturated.'''

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values at or above saturation limit
    data.data[0, 5, 5] = 60000
    data.data[1, 5, 5] = 62000
    data.data[2, 5, 5] = 62000
    data.data[3, 5, 5] = 60000
    data.data[4, 5, 5] = 62000

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Make sure all groups are flagged
    assert np.all(output.groupdq[:, 5, 5] == dqflags.group['SATURATED'])

def test_dq_propagation(setup_wfi_datamodels):
    '''Check PIXELDQ propagation.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    dqval1 = 5
    dqval2 = 10

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add DQ values to the data and reference file
    data.pixeldq[5, 5] = dqval1
    satmap.dq[5, 5] = dqval2

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Make sure DQ values from data and reference file are added in the output
    assert output.pixeldq[5, 5] == dqval1 + dqval2


def test_no_sat_check(setup_wfi_datamodels):
    '''Check that pixels flagged with NO_SAT_CHECK in the reference file get
       added to the DQ mask and are not flagged as saturated.'''

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 5, 5] = 10
    data.data[1, 5, 5] = 20000
    data.data[2, 5, 5] = 40000
    data.data[3, 5, 5] = 60000
    data.data[4, 5, 5] = 62000   # Signal reaches saturation limit

    # Set saturation value in the saturation model & DQ value for NO_SAT_CHECK
    satmap.data[5, 5] = satvalue
    satmap.dq[5, 5] = dqflags.pixel['NO_SAT_CHECK']

    # Also set an existing DQ flag in input science data
    data.pixeldq[5, 5] = dqflags.pixel['DO_NOT_USE']

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Make sure output GROUPDQ does not get flagged as saturated
    # Make sure PIXELDQ is set to NO_SAT_CHECK and original flag
    assert np.all(output.groupdq[:, 5, 5] != dqflags.group['SATURATED'])
    # Test that saturation bit is NOT set
    assert np.all(output.groupdq[:, 5, 5] & (1 << dqflags.group['SATURATED'].bit_length()-1) == 0)
    assert output.pixeldq[5, 5] == (dqflags.pixel['NO_SAT_CHECK'] +
                                    dqflags.pixel['DO_NOT_USE'])


def test_nans_in_mask(setup_wfi_datamodels):
    '''Check that pixels in the reference files that have value NaN are not
       flagged as saturated in the data and that in the PIXELDQ array the
       pixel is set to NO_SAT_CHECK.'''

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20

    data, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 5, 5] = 10
    data.data[1, 5, 5] = 20000
    data.data[2, 5, 5] = 40000
    data.data[3, 5, 5] = 60000
    data.data[4, 5, 5] = 62000

    # Set saturation value for pixel to NaN
    satmap.data[5, 5] = np.nan

    # Run the pipeline
    output = flag_saturation(data, satmap)

    # Check that output GROUPDQ is not flagged as saturated
    assert np.all(output.groupdq[:, 5, 5] != dqflags.group['SATURATED'])
    # Check that output PIXELDQ is set to NO_SAT_CHECK
    assert output.pixeldq[5, 5] == dqflags.pixel['NO_SAT_CHECK']

@pytest.fixture(scope='function')
def setup_wfi_datamodels():
    ''' Set up fake WFI data to test.'''

    def _models(ngroups, nrows, ncols):

        # Create ramp data
        data_model = testutil.mk_ramp(shape=(ngroups, nrows, ncols))

        # Create saturation reference data
        saturation_model = testutil.mk_saturation(shape=(nrows, ncols))

        return data_model, saturation_model

    return _models
