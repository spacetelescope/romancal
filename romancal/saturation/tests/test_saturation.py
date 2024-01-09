"""

Unit tests for saturation flagging

"""

import numpy as np
import pytest
from astropy import units as u
from roman_datamodels import maker_utils
from roman_datamodels.datamodels import ScienceRawModel

from romancal.lib import dqflags
from romancal.saturation import SaturationStep
from romancal.saturation.saturation import flag_saturation


def test_basic_saturation_flagging(setup_wfi_datamodels):
    """Check that the saturation flag is set when a pixel value is above the
    threshold given by the reference file."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000
    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 0 * ramp.data.unit
    ramp.data[1, 5, 5] = 20000 * ramp.data.unit
    ramp.data[2, 5, 5] = 40000 * ramp.data.unit
    ramp.data[3, 5, 5] = 60000 * ramp.data.unit  # Signal reaches saturation limit
    ramp.data[4, 5, 5] = 62000 * ramp.data.unit

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Make sure that groups with signal > saturation limit get flagged
    satindex = np.argmax(output.data.value[:, 5, 5] == satvalue)
    assert np.all(output.groupdq[satindex:, 5, 5] == dqflags.group["SATURATED"])


def test_read_pattern_saturation_flagging(setup_wfi_datamodels):
    """Check that the saturation threshold varies depending on how the reads
    are allocated into resultants."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000
    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 0 * ramp.data.unit
    ramp.data[1, 5, 5] = 20000 * ramp.data.unit
    ramp.data[2, 5, 5] = 40000 * ramp.data.unit
    ramp.data[3, 5, 5] = 60000 * ramp.data.unit  # Signal reaches saturation limit
    ramp.data[4, 5, 5] = 62000 * ramp.data.unit

    # set read_pattern to have many reads in the third resultant, so that
    # its mean exposure time is much smaller than its last read time
    # (in this case, the ratio is 13 / 20).
    # This means that the effective saturation for the third resultant
    # is 60000 * 13 / 20 = 39000 and the third resultant should be marked
    # saturated.
    ramp.meta.exposure.read_pattern = [
        [1],
        [2],
        [3, 4, 5, 6, 7, 8, 9, 10],
        [11],
        [12],
        [13],
    ]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Make sure that groups after the third get flagged
    assert np.all(output.groupdq[2:, 5, 5] == dqflags.group["SATURATED"])


def test_ad_floor_flagging(setup_wfi_datamodels):
    """Check that the ad_floor flag is set when a pixel value is zero or
    negative."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 0 * ramp.data.unit  # Signal at bottom rail - low saturation
    ramp.data[1, 5, 5] = 0 * ramp.data.unit  # Signal at bottom rail - low saturation
    ramp.data[2, 5, 5] = 20 * ramp.data.unit
    ramp.data[3, 5, 5] = 40 * ramp.data.unit
    ramp.data[4, 5, 5] = 60 * ramp.data.unit

    # frames that should be flagged as saturation (low)
    satindxs = [0, 1]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Check if the right frames are flagged as saturated
    assert np.all(
        output.groupdq[satindxs, 5, 5]
        == dqflags.group["DO_NOT_USE"] | dqflags.group["AD_FLOOR"]
    )


def test_ad_floor_and_saturation_flagging(setup_wfi_datamodels):
    """
    Check that the ad_floor flag is set when a pixel value is zero or
    negative and the saturation flag when the pixel is above the saturation threshold.
    """

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 0 * ramp.data.unit  # Signal at bottom rail - low saturation
    ramp.data[1, 5, 5] = 0 * ramp.data.unit  # Signal at bottom rail - low saturation
    ramp.data[2, 5, 5] = 20 * ramp.data.unit
    ramp.data[3, 5, 5] = 40 * ramp.data.unit
    ramp.data[4, 5, 5] = 61000 * ramp.data.unit  # Signal above the saturation threshold

    # frames that should be flagged as ad_floor
    floorindxs = [0, 1]
    # frames that should be flagged as saturation
    satindxs = [4]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Check if the right frames are flagged as ad_floor
    assert np.all(
        output.groupdq[floorindxs, 5, 5]
        == dqflags.group["DO_NOT_USE"] | dqflags.group["AD_FLOOR"]
    )
    # Check if the right frames are flagged as saturated
    assert np.all(output.groupdq[satindxs, 5, 5] == dqflags.group["SATURATED"])


def test_signal_fluctuation_flagging(setup_wfi_datamodels):
    """Check that once a pixel is flagged as saturated in a group, all
    subsequent groups should also be flagged as saturated, even if the
    signal value drops back below saturation."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 10 * ramp.data.unit
    ramp.data[1, 5, 5] = 20000 * ramp.data.unit
    ramp.data[2, 5, 5] = 40000 * ramp.data.unit
    ramp.data[3, 5, 5] = 60000 * ramp.data.unit  # Signal reaches saturation limit
    ramp.data[4, 5, 5] = 40000 * ramp.data.unit  # Signal drops below saturation limit

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Make sure that all groups after first saturated group are flagged
    satindex = np.argmax(output.data.value[:, 5, 5] == satvalue)
    assert np.all(output.groupdq[satindex:, 5, 5] == dqflags.group["SATURATED"])


def test_all_groups_saturated(setup_wfi_datamodels):
    """Check case where all groups are saturated."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values at or above saturation limit
    ramp.data[0, 5, 5] = 60000 * ramp.data.unit
    ramp.data[1, 5, 5] = 62000 * ramp.data.unit
    ramp.data[2, 5, 5] = 62000 * ramp.data.unit
    ramp.data[3, 5, 5] = 60000 * ramp.data.unit
    ramp.data[4, 5, 5] = 62000 * ramp.data.unit

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Make sure all groups are flagged
    assert np.all(output.groupdq[:, 5, 5] == dqflags.group["SATURATED"])


def test_dq_propagation(setup_wfi_datamodels):
    """Check PIXELDQ propagation."""

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    dqval1 = 5
    dqval2 = 10

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add DQ values to the data and reference file
    ramp.pixeldq[5, 5] = dqval1
    satmap.dq[5, 5] = dqval2

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Make sure DQ values from data and reference file are added in the output
    assert output.pixeldq[5, 5] == dqval1 + dqval2


def test_no_sat_check(setup_wfi_datamodels):
    """Check that pixels flagged with NO_SAT_CHECK in the reference file get
    added to the DQ mask and are not flagged as saturated."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 10 * ramp.data.unit
    ramp.data[1, 5, 5] = 20000 * ramp.data.unit
    ramp.data[2, 5, 5] = 40000 * ramp.data.unit
    ramp.data[3, 5, 5] = 60000 * ramp.data.unit
    ramp.data[4, 5, 5] = 62000 * ramp.data.unit  # Signal reaches saturation limit

    # Set saturation value in the saturation model & DQ value for NO_SAT_CHECK
    satmap.data[5, 5] = satvalue * satmap.data.unit
    satmap.dq[5, 5] = dqflags.pixel["NO_SAT_CHECK"]

    # Also set an existing DQ flag in input science data
    ramp.pixeldq[5, 5] = dqflags.pixel["DO_NOT_USE"]

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Make sure output GROUPDQ does not get flagged as saturated
    # Make sure PIXELDQ is set to NO_SAT_CHECK and original flag
    assert np.all(output.groupdq[:, 5, 5] != dqflags.group["SATURATED"])
    # Test that saturation bit is NOT set
    assert np.all(
        output.groupdq[:, 5, 5] & (1 << dqflags.group["SATURATED"].bit_length() - 1)
        == 0
    )
    assert output.pixeldq[5, 5] == (
        dqflags.pixel["NO_SAT_CHECK"] + dqflags.pixel["DO_NOT_USE"]
    )


def test_nans_in_mask(setup_wfi_datamodels):
    """Check that pixels in the reference files that have value NaN are not
    flagged as saturated in the data and that in the PIXELDQ array the
    pixel is set to NO_SAT_CHECK."""

    # Create inputs, and data and saturation models
    ngroups = 5
    nrows = 20
    ncols = 20

    ramp, satmap = setup_wfi_datamodels(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    ramp.data[0, 5, 5] = 10 * ramp.data.unit
    ramp.data[1, 5, 5] = 20000 * ramp.data.unit
    ramp.data[2, 5, 5] = 40000 * ramp.data.unit
    ramp.data[3, 5, 5] = 60000 * ramp.data.unit
    ramp.data[4, 5, 5] = 62000 * ramp.data.unit

    # Set saturation value for pixel to NaN
    satmap.data[5, 5] = np.nan * satmap.data.unit

    # Run the pipeline
    output = flag_saturation(ramp, satmap)

    # Check that output GROUPDQ is not flagged as saturated
    assert np.all(output.groupdq[:, 5, 5] != dqflags.group["SATURATED"])
    # Check that output PIXELDQ is set to NO_SAT_CHECK
    assert output.pixeldq[5, 5] == dqflags.pixel["NO_SAT_CHECK"]


def test_saturation_getbestref(setup_wfi_datamodels):
    """Check that when CRDS returns N/A for the reference file the
    step is skipped"""

    # Set test size
    shape = (2, 20, 20)

    # Create test science raw model
    wfi_sci_raw = maker_utils.mk_level1_science_raw(shape=shape)
    wfi_sci_raw.meta.instrument.name = "WFI"
    wfi_sci_raw.meta.instrument.detector = "WFI01"
    wfi_sci_raw.meta.instrument.optical_element = "F158"
    wfi_sci_raw.meta["guidestar"]["gw_window_xstart"] = 1012
    wfi_sci_raw.meta["guidestar"]["gw_window_xsize"] = 16
    wfi_sci_raw.meta.exposure.type = "WFI_IMAGE"
    wfi_sci_raw.data = u.Quantity(
        np.ones(shape, dtype=np.uint16), u.DN, dtype=np.uint16
    )
    wfi_sci_raw_model = ScienceRawModel(wfi_sci_raw, dq=True)

    # Run the pipeline
    result = SaturationStep.call(wfi_sci_raw_model, override_saturation="N/A")
    assert result.meta.cal_step.saturation == "SKIPPED"


@pytest.fixture(scope="function")
def setup_wfi_datamodels():
    """Set up fake WFI data to test."""

    def _models(ngroups, nrows, ncols):
        # Create ramp data
        ramp_model = maker_utils.mk_ramp(shape=(ngroups, nrows, ncols))

        # Create saturation reference data
        saturation_model = maker_utils.mk_saturation(shape=(nrows, ncols))

        return ramp_model, saturation_model

    return _models
