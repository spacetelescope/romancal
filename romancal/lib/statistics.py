import logging

import numpy as np
from astropy.stats import mad_std
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

    if not getattr(model.meta, "statistics", False):
        logger.debug("Creating meta.statistics node...")
        model.meta.statistics = {}

    # initialize statistics with default values
    stats = {
        "zodiacal_light": -1.0,
        "image_median": np.nan,
        "image_rms": np.nan,
        "good_pixel_fraction": np.nan,
    }

    if (
        hasattr(model, "data")
        and model.data is not None
        and not np.all(np.isnan(model.data))
    ):
        stats["image_median"] = float(np.nanmedian(model.data))
        stats["image_rms"] = mad_std(model.data, ignore_nan=True)
        if hasattr(model, "dq") and model.dq is not None:
            num_good = np.sum((model.dq & pixel.DO_NOT_USE) == 0)
            stats["good_pixel_fraction"] = num_good / model.data.size

    for key, value in stats.items():
        model.meta.statistics[key] = value
