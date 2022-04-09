#! /usr/bin/env python

from romancal.stpipe import RomanStep
from romancal.photom import photom
import roman_datamodels as rdm

__all__ = ["PhotomStep"]


class PhotomStep(RomanStep):
    """
    PhotomStep: Module for loading photometric conversion information from
        reference files and attaching to the input science data model
    """

    spec = """
        source_type = string(default=None)  # Process as specified source type.
    """
    reference_file_types = ['photom']

    def process(self, input):
        """Perform the photometric calibration step

        Parameters
        ----------
        input : Roman level 2 image datamodel (wfi_image-1.x.x)
            input roman datamodel

        Returns
        -------
        output_model : Roman level 2 image datamodel (wfi_image-1.x.x)
            output roman datamodel
        """

        # Open the input data model
        with rdm.open(input) as input_model:

            # Get reference file
            reference_file_names = {}
            reffile = self.get_reference_file(input_model, "photom")
            reference_file_names['photom'] = reffile if reffile != 'N/A' else None

            # Create storage for reference files as datamodels
            reference_file_models = {}

            # Open the relevant reference files as datamodels
            if reffile is not None:
                # If there are reference files, perform photom application
                reference_file_models['photom'] = rdm.open(reffile)
                self.log.debug(f'Using PHOTOM ref file: {reffile}')

                # Do the correction
                output_model = photom.apply_photom(input_model, **reference_file_models,)
                output_model.meta.cal_step.photom = 'COMPLETE'

            else:
                # Skip Photom step if no photom file
                self.log.warning('No PHOTOM reference file found')
                self.log.warning('Photom step will be skipped')
                output_model = input_model.copy()
                output_model.meta.cal_step.photom = 'SKIPPED'
                return output_model

        # Close the input and reference files
        input_model.close()
        try:
            for model in reference_file_models.values():
                model.close()
        except AttributeError:
            pass

        return output_model
