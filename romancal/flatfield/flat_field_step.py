"""
Flat-field a science image
"""

from crds.core.exceptions import CrdsLookupError
from ..stpipe import RomanStep
from . import flat_field
import roman_datamodels as rdm

__all__ = ["FlatFieldStep"]


class FlatFieldStep(RomanStep):
    """Flat-field a science image using a flatfield reference image.
    """

    reference_file_types = ["flat"]

    def process(self, step_input):

        input_model = rdm.open(step_input)

        # Get reference file paths
        try:
            reference_file_name = self.get_reference_file(input_model, "flat")
        except CrdsLookupError:
            reference_file_name = None

        if reference_file_name is not None:
            reference_file_model = rdm.open(reference_file_name)
            self.log.debug(f'Using FLAT ref file: {reference_file_name}')
        else:
            reference_file_model = None
            self.log.debug('Using FLAT ref file')

        # Do the flat-field correction
        output_model = flat_field.do_correction(
            input_model,
            reference_file_model,
        )

        # Close the input and reference files
        input_model.close()
        try:
            reference_file_model.close()
        except AttributeError:
            pass

        if self.save_results:
            try:
                self.suffix = 'flat'
            except AttributeError:
                self['suffix'] = 'flat'

        return output_model
