import logging

import numpy as np
from astropy.stats import mad_std
from roman_datamodels.datamodels import MosaicModel
from roman_datamodels.dqflags import pixel

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def populate_statistics(model):
    """
    Populate the meta.statistics attribute of the given model with image statistics.

    This method computes and assigns the following statistics to model.meta.statistics:
      - zodiacal_light (placeholder value)
      - image_median (median of the data array, ignoring NaNs)
      - image_rms (MAD standard deviation, ignoring NaNs)
      - good_pixel_fraction (fraction of pixels not flagged as DO_NOT_USE)

    Notes:
      - If the model lacks a data array, the method will skip statistics population;
      - If model.meta.statistics does not exist, it is created as an empty dictionary.

    Parameters
    ----------
    model : ImageModel or MosaicModel
        The data model object to populate statistics for (must have a 'model.data' attribute).

    Returns
    -------
    None
    """

    if not hasattr(model.meta, "statistics") or model.meta.statistics is None:
        logger.debug("Creating meta.statistics node...")
        model.meta.statistics = {}

    # initialize statistics with default values
    stats = {
        "zodiacal_light": -1.0,
        "image_median": np.nan,
        "image_rms": np.nan,
        "good_pixel_fraction": 0.0,
    }
    if isinstance(model, MosaicModel):
        stats.pop("zodiacal_light", None)

    if hasattr(model, "data") and model.data is not None:
        good = np.isfinite(model.data)
        if hasattr(model, "dq") and model.dq is not None:
            good &= (model.dq & pixel.DO_NOT_USE) == 0
        if hasattr(model, "err") and model.err is not None:
            good &= model.err > 0
        if not np.any(good):
            logger.warning("No good pixels found for statistics calculation.")
        else:
            stats["image_median"] = float(np.median(model.data[good]))
            stats["image_rms"] = float(mad_std(model.data[good]))
            stats["good_pixel_fraction"] = np.sum(good) / model.data.size

    for key, value in stats.items():
        model.meta.statistics[key] = value
