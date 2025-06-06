"""Apply the flux scaling"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels import datamodels

from ..datamodels import ModelLibrary
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

    def process(self, input):
        if isinstance(input, datamodels.DataModel):
            input_models = ModelLibrary([input])
            single_model = True
        elif isinstance(input, str):
            # either a single asdf filename or an association filename
            try:
                # association filename
                input_models = ModelLibrary(input)
                single_model = False
            except Exception:
                # single ASDF filename
                input_models = ModelLibrary([datamodels.open(input)])
                single_model = True
        elif isinstance(input, ModelLibrary):
            input_models = input
            single_model = False
        else:
            raise TypeError(
                "Input must be an ASN filename, a ModelLibrary, "
                "a single ASDF filename, or a single Roman DataModel."
            )

        with input_models:
            for index, model in enumerate(input_models):
                apply_flux_correction(model)
                model.meta.cal_step.flux = "COMPLETE"
                input_models.shelve(model, index)

        if single_model:
            with input_models:
                model = input_models.borrow(0)
                input_models.shelve(model, 0, modify=False)
            return model
        return input_models


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
