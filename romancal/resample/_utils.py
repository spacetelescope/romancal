import numpy as np
from roman_datamodels.dqflags import pixel

from romancal.lib.basic_utils import compute_var_rnoise


def compute_var_sky(model) -> None:
    """
    Add sky variance array to a datamodel.

    Parameters
    ----------
    model : `ImageModel`
        A datamodel to which the sky variance array will be added.

    Returns
    -------
    None
    """

    dnu = (model["dq"] & pixel.DO_NOT_USE) != 0
    median_data = np.median(model["data"][~dnu])
    ok_data = model["data"] != 0

    var_rnoise = compute_var_rnoise(model)
    var_sky = var_rnoise.copy()
    var_sky[ok_data] = (
        var_rnoise[ok_data]
        + model["var_poisson"][ok_data] / model["data"][ok_data] * median_data
    )

    return var_sky
