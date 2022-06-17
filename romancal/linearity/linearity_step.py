"""
Apply linearity correction to a science image
"""

from roman_datamodels import datamodels as rdd
import roman_datamodels as rdm
import numpy as np

from romancal.stpipe import RomanStep
from romancal.lib import dqflags
from stcal.linearity.linearity import linearity_correction

__all__ = ["LinearityStep"]


class LinearityStep(RomanStep):
    """
    LinearityStep: This step performs a correction for non-linear
    detector response, using the "classic" polynomial method.
    """

    reference_file_types = ['linearity']

    def process(self, input):

        # Open the input data model
        with rdm.open(input) as input_model:

            # Get the name of the linearity reference file to use
            self.lin_name = self.get_reference_file(input_model, 'linearity')
            self.log.info('Using Linearity reference file %s', self.lin_name)

            # Check for a valid reference file
            if self.lin_name == 'N/A':
                self.log.warning('No Linearity reference file found')
                self.log.warning('Linearity step will be skipped')
                result = input_model.copy()
                result.meta.cal_step['linearity'] = 'SKIPPED'

                return result

            lin_model = rdd.LinearityRefModel(self.lin_name)

            # copy poly coeffs from linearity model so Nan's can be updated
            lin_coeffs = lin_model.coeffs.copy()
            lin_dq = lin_model.dq  # 2D pixeldq from linearity model

            gdq = input_model.groupdq   # groupdq array of input model
            pdq = input_model.pixeldq   # pixeldq array of input model

            gdq = gdq[np.newaxis, :]

            output_model = input_model.copy()
            output_model.data = output_model.data[np.newaxis, :]

            # Call linearity correction function in stcal
            # The third return value is the procesed zero frame which
            # Roman does not use.
            new_data, new_pdq, _ = linearity_correction(output_model.data,
                                                     gdq, pdq, lin_coeffs,
                                                     lin_dq, dqflags.pixel)

            output_model.data = new_data[0, :, :, :]
            output_model.pixeldq = new_pdq

            # Close the reference file and update the step status
            lin_model.close()
            output_model.meta.cal_step['linearity'] = 'COMPLETE'

        if self.save_results:
            try:
                self.suffix = 'linearity'
            except AttributeError:
                self['suffix'] = 'linearity'
        return output_model
