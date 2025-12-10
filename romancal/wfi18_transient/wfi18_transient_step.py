#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import roman_datamodels as rdm

from romancal.stpipe import RomanStep
from romancal.wfi18_transient.wfi18_transient import correct_anomaly

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["WFI18TransientStep"]

log = logging.getLogger(__name__)


class WFI18TransientStep(RomanStep):
    """
    Correct a transient anomaly in the first read.

    This correction applies to detector WFI18 only.
    """

    class_alias = "wfi18_transient"
    reference_file_types: ClassVar = []

    spec = """
        mask_rows = boolean(default = False) # Mask the affected rows instead of fitting
    """

    def process(self, input_data):
        if isinstance(input_data, rdm.DataModel):
            input_model = input_data
        else:
            # Open the input data model
            input_model = rdm.open(input_data)

        # Ignore any data that is not WFI18
        if input_model.meta.instrument.detector != "WFI18":
            input_model.meta.cal_step.wfi18_transient = "N/A"
            return input_model

        log.info("Correcting the first read transient anomaly for WFI18")
        output_model = correct_anomaly(input_model, mask_rows=self.mask_rows)
        output_model.meta.cal_step.wfi18_transient = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "wfi18_transient"
            except AttributeError:
                self["suffix"] = "wfi18_transient"

        return output_model
