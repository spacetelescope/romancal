"""Apply the flux scaling"""
import logging

import astropy.units as u

from ..datamodels import ModelContainer
from ..stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["FluxStep"]


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

    class_alias = 'flux'

    spec = """
    """ # noqa: E501

    reference_file_types = []

    def process(self, input):
        apply_flux_correction(input)

        return input


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
    # Define the various variance arrays
    VARIANCES = ('var_rnoise', 'var_poisson', 'var_flat')

    # Check for units. Must be election/second. Otherwise, it is unknown how to
    # convert.
    if model.data.unit != u.electron / u.s:
        log.debug('Input data is not in units of e/s. Flux correction will not be done.')
        log.debug('Input data units are %s', model.data.unit)
        return

    # Apply the correction
    # Assignments into the model are done through `_instance` to avoid
    # validation errors on the units.
    log.debug('Flux correction being applied')
    c_mj = model.meta.photometry.conversion_megajanskys
    model._instance['data'] = model.data * c_mj
    for variance in VARIANCES:
        model._instance[variance] = getattr(model, variance) * c_mj**2
