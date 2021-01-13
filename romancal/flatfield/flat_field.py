#
#  Module for applying flat fielding
#

import logging
import math

import numpy as np
from gwcs.wcstools import grid_from_bounding_box

from .. import datamodels
from .. datamodels import dqflags
from .. lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

MICRONS_100 = 1.e-4  # 100 microns, in meters


def do_correction(input_model, flat=None):
    """Flat-field a Roman data model using a flat-field model

    Parameters
    ----------
    input_model : Roman data model
        Input science data model to be flat-fielded.

    flat : Roman data model, or None
        Data model containing flat-field for all instruments

    Returns
    -------
    output_model : data model
        The data model for the flat-fielded science data.
    """

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    do_flat_field(output_model, flat)

    return output_model


def do_flat_field(output_model, flat_model):
    """Apply flat-fielding, and update the output model.

    Parameters
    ----------
    output_model : Roman data model
        flat-fielded input science data model, modified in-place

    flat_model : Roman data model
        data model containing flat-field
    """

    any_updated = False  # will set True if any flats applied

    # Check to see if flat data array is smaller than science data
    if (output_model.data.shape[-1] > flat_model.data.shape[-1]) or \
       (output_model.data.shape[-2] > flat_model.data.shape[-2]):
        log.warning('Reference data array is smaller than science data')
        log.warning('Step will be skipped')

    # Apply flat to all other models
    else:
        apply_flat_field(output_model, flat_model)
        any_updated = True

    if any_updated:
        output_model.meta.cal_step.flat_field = 'COMPLETE'
    else:
        output_model.meta.cal_step.flat_field = 'SKIPPED'


def apply_flat_field(science, flat):
    """Flat field the data and error arrays.

    Extended summary
    ----------------
    The science data and error arrays will be divided by the flat field.
    The data quality array will be updated based on bad pixels in flat
    field arrays. Applies portion of flat field corresponding to science
    image subarray.

    Parameters
    ----------
    science : Roman data model
        input science data model

    flat : Roman data model
        flat field data model
    """

    # Extract subarray from reference data, if necessary
    if reffile_utils.ref_matches_sci(science, flat):
        flat_data = flat.data
        flat_dq = flat.dq
        flat_err = flat.err
    else:
        log.info("Extracting matching subarray from flat")
        sub_flat = reffile_utils.get_subarray_model(science, flat)
        flat_data = sub_flat.data.copy()
        flat_dq = sub_flat.dq.copy()
        flat_err = sub_flat.err.copy()
        sub_flat.close()

    # Find pixels in the flat that have a value of NaN and set
    # their DQ to NO_FLAT_FIELD
    flat_nan = np.isnan(flat_data)
    flat_dq[flat_nan] = np.bitwise_or(flat_dq[flat_nan],
                                      dqflags.pixel['NO_FLAT_FIELD'])

    # Find pixels in the flat that have a value of zero, and set
    # their DQ to NO_FLAT_FIELD
    flat_zero = np.where(flat_data == 0.)
    flat_dq[flat_zero] = np.bitwise_or(flat_dq[flat_zero],
                                       dqflags.pixel['NO_FLAT_FIELD'])

    # Find all pixels in the flat that have a DQ value of NO_FLAT_FIELD
    flat_bad = np.bitwise_and(flat_dq, dqflags.pixel['NO_FLAT_FIELD'])

    # Reset the flat value of all bad pixels to 1.0, so that no
    # correction is made
    flat_data[np.where(flat_bad)] = 1.0

    # Now let's apply the correction to science data and error arrays.  Rely
    # on array broadcasting to handle the cubes
    science.data /= flat_data

    # Update the variances using BASELINE algorithm.  For guider data, it has
    # not gone through ramp fitting so there is no Poisson noise or readnoise
    if not isinstance(science, datamodels.GuiderCalModel):
        flat_data_squared = flat_data**2
        science.var_poisson /= flat_data_squared
        science.var_rnoise /= flat_data_squared
        science.var_flat = science.data**2 / flat_data_squared * flat_err**2
        science.err = np.sqrt(science.var_poisson + science.var_rnoise + science.var_flat)
    else:
        flat_data_squared = flat_data**2
        science.var_flat = science.data**2 / flat_data_squared * flat_err**2
        science.err = np.sqrt(science.var_flat)

    # Combine the science and flat DQ arrays
    science.dq = np.bitwise_or(science.dq, flat_dq)
