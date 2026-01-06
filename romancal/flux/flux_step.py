"""Apply the flux scaling"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from romancal.datamodels.fileio import open_dataset

from ..stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["FluxStep"]


class FluxStep(RomanStep):
    """Apply flux scaling to count-rate data

    Parameters
    -----------
    input : str, `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.library.ModelLibrary`
        If a string is provided, it should correspond to either a single ASDF filename
        or an association filename. Alternatively, a single DataModel instance can be
        provided instead of an ASDF filename. Multiple files can be processed via
        either an association file or wrapped by a
        `~romancal.datamodels.library.ModelLibrary`.

    Returns
    -------
    output_models : `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.library.ModelLibrary`
        The models with flux applied.


    Notes
    -----
    Currently, the correction is done in-place; the inputs are directly modified if in-memory DataModels are input.
    """

    class_alias = "flux"

    spec = """
    """

    reference_file_types: ClassVar = []

    def process(self, dataset):
        input_models, dataset_type = open_dataset(
            dataset, return_type=True, as_library=True
        )
        return_lib = dataset_type in ("ModelLibrary", "asn")

        with input_models:
            for index, model in enumerate(input_models):
                apply_flux_correction(model)
                model.meta.cal_step.flux = "COMPLETE"
                input_models.shelve(model, index)

        if return_lib:
            return input_models

        with input_models:
            model = input_models.borrow(0)
            input_models.shelve(model, modify=False)
        return model


def apply_flux_correction(model):
    """Apply the flux correction

    The input model is expected to be Roman ImageModel-like.
    The photometry information is taken from `meta.photometry`.

    The model is converted in-place.

    Parameters
    ----------
    model : `ImageModel`
        The model to apply the flux correction to. The model is modified in-place.

    Notes
    -----
    The modifications to the model can result in validation issues due to change of units.
    """
    # Define the various arrays to be converted.
    DATA = ("data", "err")
    VARIANCES = ("var_rnoise", "var_poisson")
    if hasattr(model, "var_flat"):
        VARIANCES = (*VARIANCES, "var_flat")

    if model.meta.cal_step["flux"] == "COMPLETE":
        message = (
            "Input data is already in flux units of MJy/sr."
            "\nFlux correction already applied."
        )
        log.info(message)
        return

    # Apply the correction.
    # The end goal in units is to have MJy/sr. The scale is in MJy/sr also.
    # Hence the extra factor of s/DN must be applied to cancel DN/s.
    log.debug("Flux correction being applied")
    c_mj = model.meta.photometry.conversion_megajanskys
    for data in DATA:
        model[data] *= c_mj
    for variance in VARIANCES:
        model[variance] *= c_mj**2
