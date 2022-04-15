#!/usr/bin/env python
from os.path import basename
import logging
from roman_datamodels import datamodels as rdd
from ..stpipe import RomanPipeline

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

        # open the input file
        input = rdd.open(input)

        # propagate output_dir to steps that might need it
        # self.dark_current.output_dir = self.output_dir
        # self.ramp_fit.output_dir = self.output_dir

        log.debug('Exposure Processing a WFI exposure')

        self.dq_init.suffix = 'dq_init'
        result = self.dq_init(input)
        if input_filename:
            result.meta.filename = input_filename
        result = self.saturation(result)
        result = self.linearity(result)
        result = self.dark_current(result)
        result = self.jump(result)
        result = self.rampfit(result)
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
