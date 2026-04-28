#  Module for 2d saturation
#
import logging

import numpy as np
from roman_datamodels.dqflags import pixel
from stcal.saturation.saturation import flag_saturated_pixels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

ATOD_LIMIT = 65535.0  # Hard DN limit of 16-bit A-to-D converter


def flag_saturation(input_model, ref_model, n_pix_grow_sat=0, backup=0):
    """
    Flag saturated pixels, wrapping stcal.  Optionally extend flagging backwards
    in order to be conservative.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.RampModel`
        The input science data to be corrected

    ref_model : `~roman_datamodels.datamodels.SaturationRefModel`
        Saturation reference file data model

    n_pix_grow_sat : int
        Number of pixels to grow each saturated pixel by.

    backup : int
        Number of resultants to flag before each saturated resultant.

    Returns
    -------
    output_model : `~roman_datamodels.datamodels.RampModel`
        Data model with saturation, A/D floor, and do not use flags set in
        the GROUPDQ array
        The input model is modified in place and returned as the output model.
    """

    data = input_model.data[np.newaxis, :]

    # Modify input_model in place.
    gdq = input_model.groupdq[np.newaxis, :]
    pdq = input_model.pixeldq[np.newaxis, :]

    # Copy information from saturation reference file
    sat_thresh = ref_model.data.copy()
    sat_dq = ref_model.dq.copy()

    read_pattern = input_model.meta.exposure.read_pattern

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
        pixel,
        n_pix_grow_sat=n_pix_grow_sat,
        read_pattern=read_pattern,
    )

    # Save the flags in the output GROUPDQ array
    input_model.groupdq = gdq_new[0, :]

    # Save the NO_SAT_CHECK flags in the output PIXELDQ array
    input_model.pixeldq = pdq_new[0, :]

    # back saturation flagging up some frames to be safe since if the
    # non-linearity curve is sharp enough
    # the existing algorithm can fail on a large group
    # important to run this in ascending order
    for _ in range(backup):
        for i in range(len(read_pattern) - 1):
            if len(read_pattern[i]) > 1:
                input_model.groupdq[i, :, :] |= (
                    input_model.groupdq[i + 1, :, :] & pixel.SATURATED
                )

    input_model.meta.cal_step.saturation = "COMPLETE"

    return input_model
