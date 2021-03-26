""" Set DQ mask dynamically

    Parameters
    -----------
    input : model
        Filename for a FITS image or association table, or a `CubeModel`.

    output: mask
        data quality mask

"""

from .dqflags import pixel, group
from stcal import dynamicdq

def dynamic_mask(input_model, mnemonic_map=pixel):
    """
    Return a mask model given a mask with dynamic DQ flags.

    Dynamic flags define what each plane refers to using the DQ_DEF extension.

    Parameters
    ----------
    input_model : ``MaskModel``
        An instance of a Mask model defined in jwst or romancal.
    mnemonic_map : dict

    Returns
    -------
    dqmask : ndarray
        A Numpy array
    """

    return dynamicdq.dynamic_mask(input_model, mnemonic_map)

