from ..stpipe import RomanStep
from .. import datamodels
from . import flat_field
from ..datamodels import FlatModel


__all__ = ["FlatFieldStep"]


class FlatFieldStep(RomanStep):
    """Flat-field a science image using a flatfield reference image.
    """

    reference_file_types = ["flat"]

    def process(self, input):

        input_model = datamodels.open(input)

        # Get reference file paths
        reference_file_names = {}
        reffile = self.get_reference_file(input_model, "flat")
        reference_file_names['flat'] = reffile if reffile != 'N/A' else None

        # Open the relevant reference files as datamodels
        reference_file_models = {}

        if reffile is not None:
            reference_file_models['flat'] = FlatModel(reffile)
            self.log.debug(f'Using FLAT ref file: {reffile}')
        else:
            reference_file_models['flat'] = None
            self.log.debug('Using FLAT ref file')

        # Do the flat-field correction
        output_model = flat_field.do_correction(
            input_model,
            **reference_file_models,
            )

        # Close the input and reference files
        input_model.close()
        try:
            for model in reference_file_models.values():
                model.close()
        except AttributeError:
            pass

        return output_model
