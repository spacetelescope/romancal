"""Apply the flux scaling"""

from ..datamodels import ModelContainer
from ..stpipe import RomanStep


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
    Currently, the application is done in-place; the inputs are directly modified if in-memory DataModels are input.
    """  # noqa: E501

    class_alias = 'flux'

    spec = """
    """ # noqa: E501

    reference_file_types = []

    def process(self, input):
        return input
