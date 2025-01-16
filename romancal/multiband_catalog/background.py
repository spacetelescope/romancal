import numpy as np
from roman_datamodels.datamodels import ImageModel, MosaicModel

from romancal.datamodels import ModelLibrary
from romancal.source_catalog.background import RomanBackground


def subtract_background(model, box_size=1000):
    """
    Subtract the background from the input model.

    Parameters
    ----------
    model : ImageModel or MosaicModel
        The input model to subtract the background from.

    box_size : int
        The size of the box to use for the background estimation.

    Returns
    -------
    model : ImageModel or MosaicModel
        The input model with the background subtracted.
    """
    if not isinstance(model, ImageModel | MosaicModel):
        raise ValueError("The input model must be an ImageModel or MosaicModel.")

    # Subtract the background
    mask = np.isnan(model.data)
    coverage_mask = np.logical_or(np.isnan(model.err), model.err == 0)
    bkg = RomanBackground(
        model.data,
        box_size=box_size,
        mask=mask,
        coverage_mask=coverage_mask,
    )
    model.data -= bkg.background
    return model


def subtract_background_library(library, box_size=1000):
    """
    Subtract the background from all models in the input library.

    Parameters
    ----------
    library : ModelLibrary
        The input library to subtract the background from.

    box_size : int
        The size of the box to use for the background estimation.

    Returns
    -------
    library : ModelLibrary
        The input library with the background subtracted.
    """
    if not isinstance(library, ModelLibrary):
        raise TypeError("library input must be a ModelLibrary object")

    with library:
        for model in library:
            model = subtract_background(model, box_size=box_size)
            library.shelve(model)

    return library
