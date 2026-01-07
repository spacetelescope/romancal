#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels.datamodels import SaturationRefModel

from romancal.datamodels.fileio import open_dataset
from romancal.saturation import saturation
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["SaturationStep"]

log = logging.getLogger(__name__)


class SaturationStep(RomanStep):
    """
    This Step sets saturation flags.
    """

    class_alias = "saturation"

    reference_file_types: ClassVar = ["saturation"]

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        # Get the name of the saturation reference file
        self.ref_name = self.get_reference_file(input_model, "saturation")

        # Check for a valid reference file
        if self.ref_name == "N/A":
            log.warning("No SATURATION reference file found")
            log.warning("Saturation step will be skipped")
            input_model.meta.cal_step.saturation = "SKIPPED"
            return input_model

        # Open the reference file data model
        # Test for reference file
        log.info("Using SATURATION reference file: %s", self.ref_name)
        with SaturationRefModel(self.ref_name) as ref_model:
            # Perform saturation check
            saturation.flag_saturation(input_model, ref_model)

        input_model.meta.cal_step.saturation = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "saturation"
            except AttributeError:
                self["suffix"] = "saturation"

        return input_model
