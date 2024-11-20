"""
Apply linearity correction to a science image
"""

import numpy as np
from roman_datamodels import datamodels as rdd
from roman_datamodels.dqflags import pixel
from stcal.linearity.linearity import linearity_correction

from romancal.stpipe import RomanStep

__all__ = ["LinearityStep"]


class LinearityStep(RomanStep):
    """
    LinearityStep: This step performs a correction for non-linear
    detector response, using the "classic" polynomial method.
    """

    class_alias = "linearity"

    reference_file_types = ["linearity"]

    def process(self, input):
        # Open the input data model
        if isinstance(input, rdd.DataModel):
            input_model = input
        else:
            input_model = rdd.open(input)

        # Get the name of the linearity reference file to use
        self.lin_name = self.get_reference_file(input_model, "linearity")
        self.log.info("Using LINEARITY reference file: %s", self.lin_name)

        # Check for a valid reference file
        if self.lin_name == "N/A":
            self.log.warning("No LINEARITY reference file found")
            self.log.warning("Linearity step will be skipped")
            input_model.meta.cal_step["linearity"] = "SKIPPED"

            return input_model

        with rdd.LinearityRefModel(self.lin_name, memmap=False) as lin_model:
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
                input_model.data, gdq, pdq, lin_coeffs, lin_dq, pixel
            )

            input_model.data = new_data[0, :, :, :]
            input_model.pixeldq = new_pdq

        # Update the step status
        input_model.meta.cal_step["linearity"] = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "linearity"
            except AttributeError:
                self["suffix"] = "linearity"
        return input_model
