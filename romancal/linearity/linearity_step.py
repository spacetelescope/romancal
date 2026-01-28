"""
Apply linearity correction to a science image
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels import datamodels as rdd
from roman_datamodels.dqflags import pixel
from stcal.linearity.linearity import linearity_correction

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["LinearityStep"]

log = logging.getLogger(__name__)


class LinearityStep(RomanStep):
    """
    LinearityStep: This step performs a correction for non-linear
    detector response, using the "classic" polynomial method.
    """

    class_alias = "linearity"

    reference_file_types: ClassVar = ["linearity", "inverselinearity"]

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        # Get reference file names
        self.lin_name = self.get_reference_file(input_model, "linearity")
        self.ilin_name = self.get_reference_file(input_model, "inverselinearity")
        log.info("Using LINEARITY reference file: %s", self.lin_name)
        log.info("Using INVERSELINEARITY reference file: %s", self.ilin_name)

        # Check for valid reference files
        if self.lin_name == "N/A" or self.ilin_name == "N/A":
            log.warning("No LINEARITY or INVERSELINEARITY reference file found")
            log.warning("Linearity step will be skipped")
            input_model.meta.cal_step["linearity"] = "SKIPPED"
            return input_model

        with (
            rdd.LinearityRefModel(self.lin_name, memmap=False) as lin_model,
            rdd.InverselinearityRefModel(self.ilin_name, memmap=False) as ilin_model,
        ):
            lin_coeffs = lin_model.coeffs
            lin_dq = lin_model.dq
            ilin_coeffs = ilin_model.coeffs
            read_pattern = input_model.meta.exposure.read_pattern

            gdq = input_model.groupdq[np.newaxis, :]
            pdq = input_model.pixeldq
            input_model.data = input_model.data[np.newaxis, :]

            # Call linearity correction function in stcal
            # The third return value is the processed zero frame which
            # Roman does not use.
            new_data, new_pdq, _ = linearity_correction(
                input_model.data,
                gdq,
                pdq,
                lin_coeffs,
                lin_dq,
                pixel,
                ilin_coeffs=ilin_coeffs,
                read_pattern=read_pattern,
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
