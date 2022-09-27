#! /usr/bin/env python

import logging
import roman_datamodels as rdm
from roman_datamodels.datamodels import SaturationRefModel
from romancal.stpipe import RomanStep
from romancal.saturation import saturation
from romancal.lib.basic_utils import is_fully_saturated

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["SaturationStep"]


class SaturationStep(RomanStep):
    """
    This Step sets saturation flags.
    """

    reference_file_types = ['saturation']

    def process(self, input):

        # Open the input data model
        with rdm.open(input) as input_model:

            # Get the name of the saturation reference file
            self.ref_name = self.get_reference_file(input_model, 'saturation')
            self.log.info('Using SATURATION reference file %s', self.ref_name)

            # Check for a valid reference file
            if self.ref_name == 'N/A':
                self.log.warning('No SATURATION reference file found')
                self.log.warning('Saturation step will be skipped')
                result = input_model.copy()
                result.meta.cal_step.saturation = 'SKIPPED'
                return result

            # Open the reference file data model
            ref_model = SaturationRefModel(self.ref_name)

            # Perform saturation check
            sat = saturation.flag_saturation(input_model, ref_model)

            # Close the reference file and update the step status
            ref_model.close()
            sat.meta.cal_step.saturation = 'COMPLETE'

            if self.save_results:
                try:
                    self.suffix = 'saturation'
                except AttributeError:
                    self['suffix'] = 'saturation'

        # Test for fully saturated data
        if is_fully_saturated(sat):
            log.info('All pixels are saturated. Returning a zeroed-out image.')
            # Set all subsequent steps to skipped but ramp fit to create l2 image
            for step_str in ['linearity', 'dark', 'jump', 'assign_wcs',
                             'flat_field', 'photom']:
                sat.meta.cal_step[step_str] = 'SKIPPED'

        return sat
