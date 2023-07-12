"""
Apply linearity correction to a science image
"""

import numpy as np
from astropy import units as u
from roman_datamodels import datamodels as rdd
from stcal.linearity.linearity import linearity_correction

from romancal.lib import dqflags
from romancal.stpipe import RomanStep

__all__ = ["LinearityStep"]


class LinearityStep(RomanStep):
    """
    LinearityStep: This step performs a correction for non-linear
    detector response, using the "classic" polynomial method.
    """

    reference_file_types = ["linearity"]

    def process(self, input):
        # Open the input data model
        with rdd.open(input, lazy_load=False) as input_model:
            # Get the name of the linearity reference file to use
            self.lin_name = self.get_reference_file(input_model, "linearity")
            self.log.info("Using LINEARITY reference file: %s", self.lin_name)

            # Check for a valid reference file
            if self.lin_name == "N/A":
                self.log.warning("No LINEARITY reference file found")
                self.log.warning("Linearity step will be skipped")
                input_model.meta.cal_step["linearity"] = "SKIPPED"

                return input_model

            lin_model = rdd.LinearityRefModel(self.lin_name, copy_arrays=True)

            # copy poly coeffs from linearity model so Nan's can be updated
            lin_coeffs = lin_model.coeffs
            lin_dq = lin_model.dq  # 2D pixeldq from linearity model

            gdq = input_model.groupdq  # groupdq array of input model
            pdq = input_model.pixeldq  # pixeldq array of input model

            gdq = gdq[np.newaxis, :]

            input_model.data = input_model.data[np.newaxis, :]

            # Call linearity correction function in stcal
            # The third return value is the procesed zero frame which
            # Roman does not use.
            new_data, new_pdq, _ = linearity_correction(
                input_model.data.value, gdq, pdq, lin_coeffs, lin_dq, dqflags.pixel
            )

            input_model.data = u.Quantity(
                new_data[0, :, :, :], u.DN, dtype=new_data.dtype
            )
            input_model.pixeldq = new_pdq

            # Close the reference file and update the step status
            lin_model.close()
            input_model.meta.cal_step["linearity"] = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "linearity"
            except AttributeError:
                self["suffix"] = "linearity"
        return input_model
