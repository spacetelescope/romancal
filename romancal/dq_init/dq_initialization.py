import logging
import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Guide star mode exposure types
GUIDER_LIST = ['WFI_WIM_ACQ','WFI_WIM_TRACK','WFI_WSM_ACQ1','WFI_WSM_ACQ2','WFI_WSM_TRACK']


def do_dqinit(input_model, mask=None):
    """Check that the input model pixeldq attribute has the same dimensions as
    the image plane of the input model science data, call apply_dqinit, and update
    log and cal_step.

    Parameters
    ----------
    input_model : input Roman datamodel
        The Roman datamodel to be data quality  corrected

    mask_model : Roman data model, or None
        The mask model to use in the data quality correction

    Returns
    -------
    output_model : Roman datamodel
        The data quality corrected Roman datamodel
    """

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    # Determine if mask is shapewise compatable with the input
    skip_step = False
    if (input_model.meta.exposure.type in GUIDER_LIST):
        # Check to see if the shape of the mask data array
        # is not equal to the shape of the science data
        if input_model.dq.shape != mask.dq.shape:
            skip_step = True
    else:
        # Check to see if the shape of the mask data array
        # is not equal to the shape of the science data
        if input_model.pixeldq.shape != mask.dq.shape:
            skip_step = True


    # Apply or skip dq correction
    if skip_step:
        log.warning('Mask data array is not the same '
                    'shape as the science data')
        log.warning('Step will be skipped')
        output_model.meta.cal_step.dq_init = 'SKIPPED'
    else:
        output_model = apply_dqinit(input_model, mask)
        output_model.meta.cal_step.dq_init = 'COMPLETE'

    return output_model

def apply_dqinit(science, mask):
    """Apply data quality mask to the dq or pixeldq arrays.

    Extended summary
    ----------------
    The dq attribute of the MASK model is bitwise OR'd with the
    pixeldq (or dq) attribute of the input model.

    Parameters
    ----------
    science : Roman data model
        input science data model

    mask : Roman data model
        data quality mask data model
    """

    # Set model-specific data quality in output
    if (science.meta.exposure.type in GUIDER_LIST):
        science.dq = np.bitwise_or(science.dq, mask.dq)
    else:
        science.pixeldq = np.bitwise_or(science.pixeldq, mask.dq)

    return science
