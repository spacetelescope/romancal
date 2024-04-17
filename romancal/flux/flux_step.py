"""Apply the flux scaling"""

import logging

import astropy.units as u
from roman_datamodels import datamodels

from ..datamodels import ModelContainer
from ..stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["FluxStep"]


# Define expected Level 2 units
LV2_UNITS = u.DN / u.s


class FluxStep(RomanStep):
    """Apply flux scaling to count-rate data

    Parameters
    -----------
    input : str, `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.container.ModelContainer`
        If a string is provided, it should correspond to either a single ASDF filename
        or an association filename. Alternatively, a single DataModel instance can be
        provided instead of an ASDF filename. Multiple files can be processed via
        either an association file or wrapped by a
        `~romancal.datamodels.container.ModelContainer`.

    Returns
    -------
    output_models : `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.container.ModelContainer`
        The models with flux applied.


    Notes
    -----
    Currently, the correction is done in-place; the inputs are directly modified if in-memory DataModels are input.
    """  # noqa: E501

    spec = """
    """  # noqa: E501

    reference_file_types = []

    def process(self, input):
        if isinstance(input, datamodels.DataModel):
            input_models = ModelContainer([input])
            single_model = True
        elif isinstance(input, str):
            # either a single asdf filename or an association filename
            try:
                # association filename
                input_models = ModelContainer(input)
                single_model = False
            except Exception:
                # single ASDF filename
                input_models = ModelContainer([input])
                single_model = True
        elif isinstance(input, ModelContainer):
            input_models = input
            single_model = False
        else:
            raise TypeError(
                "Input must be an ASN filename, a ModelContainer, "
                "a single ASDF filename, or a single Roman DataModel."
            )

        for model in input_models:
            apply_flux_correction(model)
            model.meta.cal_step.flux = "COMPLETE"

        if single_model:
            return input_models[0]
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
    VARIANCES = ("var_rnoise", "var_poisson", "var_flat")

    if model.data.unit == model.meta.photometry.conversion_megajanskys.unit:
        message = (
            f"Input data is already in flux units of {model.meta.photometry.conversion_megajanskys.unit}."
            "\nFlux correction already applied."
        )
        log.info(message)
        return

    if model.data.unit != LV2_UNITS:
        message = (
            f"Input data units {model.data.unit} are not in the expected units of {LV2_UNITS}"
            "\nAborting flux correction"
        )
        log.error(message)
        raise ValueError(message)

    # Apply the correction.
    # The end goal in units is to have MJy/sr. The scale is in MJy/sr also.
    # Hence the extra factor of s/DN must be applied to cancel DN/s.
    log.debug("Flux correction being applied")
    c_mj = model.meta.photometry.conversion_megajanskys / model.data.unit
    for data in DATA:
        model[data] = model[data] * c_mj
    for variance in VARIANCES:
        model[variance] = model[variance] * c_mj**2
