import numpy as np
from roman_datamodels.dqflags import pixel

from romancal.lib.basic_utils import compute_var_rnoise


def compute_var_sky(model):
    """
    Compute the sky variance of a datamodel.

    The sky variance is the read noise variance plus the Poisson
    variance scaled to the median data level, so it does not depend on
    the source flux in each pixel.

    Parameters
    ----------
    model : `ImageModel`
        The datamodel.

    Returns
    -------
    var_sky : np.ndarray
        The sky variance.
    """

    dnu = (model["dq"] & pixel.DO_NOT_USE) != 0
    data = np.asarray(model["data"], dtype=np.float32)
    var_poisson = np.asarray(model["var_poisson"], dtype=np.float32)
    median_data = np.nanmedian(data[~dnu])

    var_rnoise = compute_var_rnoise(model)
    if not median_data > 0:
        return var_rnoise.copy()

    # Poisson variance per unit flux. It is unreliable where the data
    # are near zero or negative (var_poisson is clipped at zero there),
    # so use its median value in those pixels.
    bright = data > 0.1 * median_data
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(bright, var_poisson / data, np.nan)
    ratio = np.where(bright, ratio, np.nanmedian(ratio[bright & ~dnu]))
    var_sky = var_rnoise + ratio * median_data

    return var_sky
