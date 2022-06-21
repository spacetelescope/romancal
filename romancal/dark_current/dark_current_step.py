#! /usr/bin/env python

import numpy as np
from romancal.stpipe import RomanStep
from roman_datamodels import datamodels as rdd
from stcal.dark_current import dark_sub
from roman_datamodels.testing import utils as testutil


__all__ = ["DarkCurrentStep"]


class DarkCurrentStep(RomanStep):
    """
    DarkCurrentStep: Performs dark current correction by subtracting
    dark current reference data from the input science data model.
    """

    spec = """
        dark_output = output_file(default = None) # Dark corrected model
    """

    reference_file_types = ['dark']

    def process(self, input):

        # Open the input data model
        with rdd.open(input) as input_model:

            # Get the name of the dark reference file to use
            self.dark_name = self.get_reference_file(input_model, 'dark')
            self.log.info('Using DARK reference file %s', self.dark_name)

            # Open dark model
            dark_model = rdd.open(self.dark_name)

            # Temporary patch to utilize stcal dark step until MA table support is fully implemented
            if 'ngroups' not in dark_model.meta.exposure:
                dark_model.meta.exposure['ngroups'] = dark_model.data.shape[0]
            if 'nframes' not in dark_model.meta.exposure:
                dark_model.meta.exposure['nframes'] = input_model.meta.exposure.nframes
            if 'groupgap' not in dark_model.meta.exposure:
                dark_model.meta.exposure['groupgap'] = input_model.meta.exposure.groupgap

            # Reshaping data variables for stcal compatibility
            input_model.data = input_model.data.astype(np.float32)[np.newaxis, :]
            input_model.groupdq = input_model.groupdq[np.newaxis, :]
            input_model.err = input_model.err[np.newaxis, :]

            # Do the dark correction
            out_data, dark_data = dark_sub.do_correction(
                input_model, dark_model, self.dark_output
            )

            # Save dark data to file
            if dark_data is not None and dark_data.save:
                save_dark_data_as_dark_model(dark_data, dark_model)
            dark_model.close()

            # Reshaping data variables back from stcal
            input_model.data = input_model.data[0]
            input_model.groupdq = input_model.groupdq[0]
            input_model.err = input_model.err[0]

            # Convert data to RampModel
            out_ramp = dark_output_data_as_ramp_model(out_data, input_model)

        if self.save_results:
            try:
                self.suffix = 'darkcurrent'
            except AttributeError:
                self['suffix'] = 'darkcurrent'
            dark_model.close()

        return out_ramp


def save_dark_data_as_dark_model(dark_data, dark_model):
    """
    Save dark data from the dark current step as the appropriate dark model.

    Parameters
    ----------
    dark_data: DarkData
        Dark data used in the dark current step.

    dark_model: DarkRefModel
        The input dark model from reference.
    """

    # Create DarkRef object and copy dark data to it
    out_dark = testutil.mk_dark(shape=dark_data.data.shape)
    out_dark.data = dark_data.data
    out_dark.dq = dark_data.groupdq
    out_dark.err = dark_data.err

    # Temporary patch to utilize stcal dark step until MA table support is fully implemented
    out_dark.meta.exposure['nframes'] = dark_data.exp_nframes
    out_dark.meta.exposure['ngroups'] = dark_data.exp_ngroups
    out_dark.meta.exposure['groupgap'] = dark_data.exp_groupgap

    # Create DarkRefModel and write to file
    out_dark_model = rdd.DarkRefModel(out_dark)
    out_dark_model.save(dark_data.output_name)
    out_dark_model.close()


def dark_output_data_as_ramp_model(out_data, input_model):
    """
    Convert computed output data from the dark step to a RampModel.

    Parameters
    ----------
    out_data: ScienceData
        Computed science data from the dark current step.

    input_model: RampModel
        The input ramp model from which to subtract the dark current.

    Return
    ------
    out_model: RampModel
        The output ramp model from the dark current step.
    """

    # Copy input model as a base to preserve everything in addition to the dark output
    out_model = input_model.copy()
    out_model.meta.cal_step.dark = out_data.cal_step

    if out_data.cal_step == "SKIPPED":
        return out_model

    # Removing integration dimension from variables (added for stcal
    # compatibility)
    # Roman 3D
    out_model.data = out_data.data[0]
    out_model.groupdq = out_data.groupdq[0]
    # Roman 2D
    out_model.pixeldq = out_data.pixeldq
    out_model.err = out_data.err[0]

    return out_model
