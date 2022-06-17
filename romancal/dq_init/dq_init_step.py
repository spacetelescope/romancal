#! /usr/bin/env python

import numpy as np

from romancal.stpipe import RomanStep
from romancal.dq_init import dq_initialization
from roman_datamodels.datamodels import RampModel
import roman_datamodels as rdm
from roman_datamodels.testing import utils as testutil

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

        # Convert to RampModel if needed
        if not (isinstance(input_model, RampModel)):
            # Create base ramp node with dummy values (for validation)
            input_ramp = testutil.mk_ramp(input_model.shape)

            # Copy input_model contents into RampModel
            for key in input_model.keys():
                # If a dictionary (like meta), overwrite entires (but keep
                # required dummy entries that may not be in input_model)
                if isinstance(input_ramp[key], dict):
                    input_ramp[key].update(input_model.__getattr__(key))
                elif isinstance(input_ramp[key], np.ndarray):
                    # Cast input ndarray as RampModel dtype
                    input_ramp[key] = input_model.__getattr__(key).astype(input_ramp[key].dtype)
                else:
                    input_ramp[key] = input_model.__getattr__(key)

            # Create model from node
            init_model = RampModel(input_ramp)
        else:
            init_model = input_model

        # Get reference file path
        reference_file_name = self.get_reference_file(init_model, "mask")

        # Test for reference file
        if reference_file_name is not None:
            # If there are mask files, perform dq step
            # Open the relevant reference files as datamodels
            reference_file_model = rdm.open(reference_file_name)
            self.log.debug(f'Using MASK ref file: {reference_file_name}')

            # Apply the DQ step
            output_model = dq_initialization.do_dqinit(
                init_model,
                reference_file_model,
            )

            # copy original border reference file arrays (data and dq)
            # to their own attributes. they will also remain attached to
            # the science data until they are trimmed at ramp_fit
            # these arrays include the overlap regions in the corners

            output_model.border_ref_pix_right = output_model.data[:, :, -4:]
            output_model.border_ref_pix_left = output_model.data[:, :, :4]
            output_model.border_ref_pix_top = output_model.data[:, :4, :]
            output_model.border_ref_pix_bottom = output_model.data[:, -4:, :]

            output_model.dq_border_ref_pix_right = output_model.pixeldq[:, -4:]
            output_model.dq_border_ref_pix_left = output_model.pixeldq[:, :4]
            output_model.dq_border_ref_pix_top = output_model.pixeldq[:4, :]
            output_model.dq_border_ref_pix_bottom = output_model.pixeldq[-4:, :]

        else:
            # Skip DQ step if no mask files
            reference_file_model = None
            self.log.warning('No MASK reference file found.')
            self.log.warning('DQ initialization step will be skipped.')

            output_model = init_model
            output_model.meta.cal_step.dq_init = 'SKIPPED'


        # Close the input and reference files
        input_model.close()
        init_model.close()
        try:
            reference_file_model.close()
        except AttributeError:
            pass

        if self.save_results:
            try:
                self.suffix = 'dqinit'
            except AttributeError:
                self['suffix'] = 'dqinit'

        return output_model
