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
        The data model object to populate statistics for (must have a 'model.meta.data' attribute).

    Returns
    -------
    None
    """
    if not hasattr(model, "data") or model.data is None:
        logger.debug("Model has no data array; skipping statistics population.")
        return

    if getattr(model.meta, "statistics", None) is None:
        logger.debug("Creating meta.statistics node...")
        model.meta.statistics = {}

    img_median = float(np.nanmedian(model.data))

    img_rms = mad_std(model.data, ignore_nan=True)

    good_pix_frac = 0.0
    if hasattr(model, "dq") and model.dq is not None:
        num_good = np.sum((model.dq & pixel.DO_NOT_USE) == 0)
        good_pix_frac = num_good / model.data.size

        model.meta.statistics.zodiacal_light = -1.0  # placeholder
        model.meta.statistics.image_median = img_median
        model.meta.statistics.image_rms = img_rms
        model.meta.statistics.good_pixel_fraction = good_pix_frac
