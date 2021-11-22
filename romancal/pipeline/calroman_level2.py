#!/usr/bin/env python
from os.path import basename
import logging
from ..stpipe import RomanPipeline
from roman_datamodels import datamodels as rdd

# step imports
from romancal.assign_wcs import AssignWcsStep
from romancal.dq_init import dq_init_step
from romancal.flatfield import FlatFieldStep
from romancal.jump import jump_step
from romancal.dark_current import DarkCurrentStep
from romancal.linearity import LinearityStep
from romancal.ramp_fitting import ramp_fit_step
from romancal.saturation import SaturationStep

__all__ = ['Level2Pipeline']

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class Level2Pipeline(RomanPipeline):
    """
    Level21Pipeline: Apply all calibration steps to raw Roman WFI
    ramps to produce a 2-D slope product. Included steps are:
    dq_init, saturation, linearity, dark currnet, jump detection, and ramp_fit.
    """

    class_alias = "calroman_level2"

    spec = """
        save_calibrated_ramp = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {'dq_init': dq_init_step.DQInitStep,
                 'saturation': SaturationStep,
                 'linearity': LinearityStep,
                 'dark_current': DarkCurrentStep,
                 'jump': jump_step.JumpStep,
                 'rampfit': ramp_fit_step.RampFitStep,
                 'assign_wcs': AssignWcsStep,
                 'flatfield': FlatFieldStep
                 }

    # start the actual processing
    def process(self, input):

        """ Process the Roman WFI data """

        log.info('Starting calroman_level2 ...')
        if isinstance(input, str):
            input_filename = basename(input)

        # open the input file
        input = rdd.open(input)

        # propagate output_dir to steps that might need it
        # self.dark_current.output_dir = self.output_dir
        # self.ramp_fit.output_dir = self.output_dir

        log.debug('Processing a WFI exposure')

        result = self.dq_init(input)
        if input_filename:
            result.meta.filename = input_filename
        result = self.saturation(result)
        result = self.linearity(result)
        result = self.dark_current(result)
        result = self.jump(result)
        result = self.rampfit(result)
        result = self.assign_wcs(result)
        result = self.flatfield(result)

        # setup output_file for saving
        self.setup_output(result)
        log.info('calroman_level2 pipeline ending...')

        return result

    def setup_output(self, input):
        """ Determine the proper file name suffix to use later """
        if input.meta.cal_step.ramp_fit == 'COMPLETE':
            self.suffix = 'rate'
            input.meta.filename = input.meta.filename.replace(
                '_uncal', '')
            input['output_file'] = input.meta.filename
            self.output_file = input.meta.filename
        else:
            self.suffix = 'ramp'
