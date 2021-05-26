import logging
import numpy as np
from stcal.dq_init import dq_initialization

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Guide star mode exposure types
GUIDER_LIST = ['WFI_WIM_ACQ','WFI_WIM_TRACK','WFI_WSM_ACQ1','WFI_WSM_ACQ2','WFI_WSM_TRACK']

def do_dqinit(input_model, mask_model):
    """Perform the dq_init step on a Roman datamodel

    Parameters
    ----------
    input_model : input Roman datamodel
        The Roman datamodel to be corrected

    mask_model : mask datamodel
        The mask model to use in the correction

    Returns
    -------
    output_model : Roman datamodel
        The corrected Roman datamodel
    """

    return dq_initialization.do_dqinit(input_model, mask_model, GUIDER_LIST)


def check_dimensions(input_model):
    """Check that the input model pixeldq attribute has the same dimensions as
    the image plane of the input model science data
    If it has dimensions (0,0), create an array of zeros with the same shape
    as the image plane of the input model.

    For the guiding modes, the GuiderRawModel has only a regular dq array (no pixeldq or groupdq)

    Parameters
    ----------
    input_model : Raw Roman datamodel
        input Roman datamodel

    Returns
    -------
    None
    """

    return dq_initialization.check_dimensions(input_model)
