"""
Module for applying flat fielding
"""

import logging
import numpy as np

from romancal.lib import dqflags

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
    if flat_model is not None and \
            output_model.data.shape != flat_model.data.shape:
        # Check to see if flat data array is smaller than science data
        log.warning('Flat data array is not the same '
                    'shape as the science data')
        log.warning('Step will be skipped')
        output_model.meta.cal_step.flat_field = 'SKIPPED'
    elif output_model.meta.exposure.type != "WFI_IMAGE":
        # Check to see if attempt to flatten non-Image data
        log.info('Skipping flat field for spectral exposure.')
        output_model.meta.cal_step.flat_field = 'SKIPPED'
    elif flat_model is None:
        # Check to see if attempt to flatten non-Image data
        log.info('Skipping flat field - no flat reference file.')
        output_model.meta.cal_step.flat_field = 'SKIPPED'
    else:
        apply_flat_field(output_model, flat_model)
        output_model.meta.cal_step.flat_field = 'COMPLETE'


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
    flat_data = flat.data.copy()
    flat_dq = flat.dq.copy()
    flat_err = flat.err.copy()
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
    flat_data_squared = flat_data**2
    science.var_poisson /= flat_data_squared
    science.var_rnoise /= flat_data_squared
    try:
        science.var_flat = science.data**2 / flat_data_squared * flat_err**2
    except AttributeError:
        science['var_flat'] = np.zeros(shape=science.data.shape,
                                       dtype=np.float32)
        science.var_flat = science.data**2 / flat_data_squared * flat_err**2
    science.err = np.sqrt(science.var_poisson +
                          science.var_rnoise + science.var_flat)

    # Combine the science and flat DQ arrays
    science.dq = np.bitwise_or(science.dq, flat_dq)
