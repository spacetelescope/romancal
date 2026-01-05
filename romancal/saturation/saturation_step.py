#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import roman_datamodels as rdm
from roman_datamodels.datamodels import SaturationRefModel

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

    _input_class = rdm.datamodels.RampModel

    class_alias = "saturation"

    reference_file_types: ClassVar = ["saturation"]

    def process(self, input):
        if isinstance(input, rdm.DataModel):
            input_model = input
        else:
            # Open the input data model
            input_model = rdm.open(input)

        # Get the name of the saturation reference file
        self.ref_name = self.get_reference_file(input_model, "saturation")

        # Check for a valid reference file
        if self.ref_name == "N/A":
            log.warning("No SATURATION reference file found")
            log.warning("Saturation step will be skipped")
            result = input_model.copy()
            result.meta.cal_step.saturation = "SKIPPED"
            return result

        # Open the reference file data model
        # Test for reference file
        log.info("Using SATURATION reference file: %s", self.ref_name)
        ref_model = SaturationRefModel(self.ref_name)

        # Perform saturation check
        sat = saturation.flag_saturation(input_model, ref_model)

        # Close the reference file and update the step status
        ref_model.close()
        sat.meta.cal_step.saturation = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "saturation"
            except AttributeError:
                self["suffix"] = "saturation"

        return sat
