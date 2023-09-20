#  Module for 2d saturation
#
import logging

import numpy as np
from stcal.saturation.saturation import flag_saturated_pixels

from romancal.lib import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ATOD_LIMIT = 65535.0  # Hard DN limit of 16-bit A-to-D converter


def flag_saturation(input_model, ref_model):
    """
    Short Summary
    -------------
    Function to call stcal for flagging saturated pixels.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        The input science data to be corrected

    ref_model : `~roman_datamodels.datamodels.SaturationRefModel`
        Saturation reference file data model

    Returns
    -------
    output_model : `~roman_datamodels.datamodels.RampModel`
        Data model with saturation, A/D floor, and do not use flags set in
        the GROUPDQ array
        The input model is modified in place and returned as the output model.
    """

    data = input_model.data[np.newaxis, :].value

    # Modify input_model in place.
    gdq = input_model.groupdq[np.newaxis, :]
    pdq = input_model.pixeldq[np.newaxis, :]

    # Copy information from saturation reference file
    sat_thresh = ref_model.data.value.copy()
    sat_dq = ref_model.dq.copy()

    # Obtain dq arrays updated for saturation
    # The third variable is the processed ZEROFRAME, which is not
    # used in romancal, so is always None.
    gdq_new, pdq_new, _ = flag_saturated_pixels(
        data,
        gdq,
        pdq,
        sat_thresh,
        sat_dq,
        ATOD_LIMIT,
        dqflags.pixel,
        n_pix_grow_sat=0,
        read_pattern=input_model.meta.exposure.read_pattern,
    )

    # Save the flags in the output GROUPDQ array
    input_model.groupdq = gdq_new[0, :]

    # Save the NO_SAT_CHECK flags in the output PIXELDQ array
    input_model.pixeldq = pdq_new[0, :]

    return input_model
