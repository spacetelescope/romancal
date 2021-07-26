#! /usr/bin/env python

from romancal.stpipe import RomanStep
from romancal.dq_init import dq_initialization
import roman_datamodels as rdm


__all__ = ["DQInitStep"]


class DQInitStep(RomanStep):
    """Initialize the Data Quality extension from the
    mask reference file.

    The dq_init step initializes the pixeldq attribute of the
    input datamodel using the MASK reference file.  For some
    Guiding and Image model types, initalize the dq attribute of
    the input model instead.  The dq attribute of the MASK model
    is bitwise OR'd with the pixeldq (or dq) attribute of the
    input model.
    """


    reference_file_types = ['mask']

    def process(self, input):
        """Perform the dq_init calibration step

        Parameters
        ----------
        input : Roman datamodel
            input roman datamodel

        Returns
        -------
        output_model : Roman datamodel
            result roman datamodel
        """
        # Open datamodel
        input_model = self.open_model(input)

        # Get reference file paths
        reference_file_names = {}
        reffile = self.get_reference_file(input_model, "mask")
        reference_file_names['mask'] = reffile if reffile != 'N/A' else None

        # Open the relevant reference files as datamodels
        reference_file_models = {}

        if reffile is not None:
            reference_file_models['mask'] = rdm.open(reffile)
            self.log.debug(f'Using MASK ref file: {reffile}')
        else:
            reference_file_models['mask'] = None
            self.log.debug('Using MASK ref file')

        # Apply the DQ step
        output_model = dq_initialization.do_dqinit(
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
