#!/usr/bin/env python
from os.path import basename
import logging

import numpy as np
from roman_datamodels import datamodels as rdd
from ..stpipe import RomanPipeline
from romancal.lib.basic_utils import is_fully_saturated
from romancal.lib import dqflags

# step imports
from romancal.assign_wcs import AssignWcsStep
from romancal.dq_init import dq_init_step
from romancal.flatfield import FlatFieldStep
from romancal.jump import jump_step
from romancal.dark_current import DarkCurrentStep
from romancal.linearity import LinearityStep
from romancal.photom import PhotomStep
from romancal.ramp_fitting import ramp_fit_step
from romancal.saturation import SaturationStep

__all__ = ['ExposurePipeline']

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)



class ExposurePipeline(RomanPipeline):
    """
    ExposurePipeline: Apply all calibration steps to raw Roman WFI
    ramps to produce a 2-D slope product. Included steps are:
    dq_init, saturation, linearity, dark current, jump detection, ramp_fit,
    and assign_wcs. The flat field step is only applied to WFI imaging data.
    """

    class_alias = "roman_elp"

    spec = """
        save_calibrated_ramp = boolean(default=False)
        save_results = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {'dq_init': dq_init_step.DQInitStep,
                 'saturation': SaturationStep,
                 'linearity': LinearityStep,
                 'dark_current': DarkCurrentStep,
                 'jump': jump_step.JumpStep,
                 'rampfit': ramp_fit_step.RampFitStep,
                 'assign_wcs': AssignWcsStep,
                 'flatfield': FlatFieldStep,
                 'photom': PhotomStep,
                 }

    # start the actual processing
    def process(self, input):

        """ Process the Roman WFI data """

        log.info('Starting Roman exposure calibration pipeline ...')
        if isinstance(input, str):
            input_filename = basename(input)
        else:
            input_filename = None

        # open the input file
        input = rdd.open(input)

        log.debug('Exposure Processing a WFI exposure')

        self.dq_init.suffix = 'dq_init'
        result = self.dq_init(input)
        if input_filename:
            result.meta.filename = input_filename
        result = self.saturation(result)

        # Test for fully saturated data
        if is_fully_saturated(result):
            log.info('All pixels are saturated. Returning a zeroed-out image.')

            # Return zeroed-out image file (stopping pipeline)
            return self.create_fully_saturated_zeroed_image(result)

        result = self.linearity(result)
        result = self.dark_current(result)
        result = self.jump(result)
        result = self.rampfit(result)

        # Test for fully saturated data
        if "groupdq" in result.keys():
            if is_fully_saturated(result):
                # Set all subsequent steps to skipped
                for step_str in ['assign_wcs', 'flat_field', 'photom']:
                    result.meta.cal_step[step_str] = 'SKIPPED'

                # Set suffix for proper output naming
                self.suffix = 'cal'

                # Return fully saturated image file (stopping pipeline)
                return result

        result = self.assign_wcs(result)
        if result.meta.exposure.type == 'WFI_IMAGE':
            result = self.flatfield(result)
        else:
            log.info('Flat Field step is being SKIPPED')
            result.meta.cal_step.flat_field = 'SKIPPED'
        result = self.photom(result)

        # setup output_file for saving
        self.setup_output(result)
        log.info('Roman exposure calibration pipeline ending...')

        return result

    def setup_output(self, input):
        """ Determine the proper file name suffix to use later """
        if input.meta.cal_step.ramp_fit == 'COMPLETE':
            self.suffix = 'cal'
            input.meta.filename = input.meta.filename.replace(
                '_uncal', '')
            input['output_file'] = input.meta.filename
            self.output_file = input.meta.filename
        else:
            self.suffix = 'ramp'

    def create_fully_saturated_zeroed_image(self, input_model):
        """
        Create zeroed-out image file
        """
        # The set order is: data, dq, var_poisson, var_rnoise, err
        fully_saturated_model = ramp_fit_step.create_image_model(input_model,
                                              (np.zeros(input_model.data.shape[1:],
                                                        dtype=input_model.data.dtype),
                                               input_model.pixeldq | input_model.groupdq[0] |
                                                                     dqflags.group['SATURATED'],
                                               np.zeros(input_model.err.shape[1:],
                                                        dtype=input_model.err.dtype),
                                               np.zeros(input_model.err.shape[1:],
                                                        dtype=input_model.err.dtype),
                                               np.zeros(input_model.err.shape[1:],
                                                        dtype=input_model.err.dtype)))

        # Set all subsequent steps to skipped
        for step_str in ['linearity', 'dark', 'jump', 'ramp_fit', 'assign_wcs',
                         'flat_field', 'photom']:
            fully_saturated_model.meta.cal_step[step_str] = 'SKIPPED'

        # Set suffix for proper output naming
        self.suffix = 'cal'

        # Return zeroed-out image file
        return fully_saturated_model
