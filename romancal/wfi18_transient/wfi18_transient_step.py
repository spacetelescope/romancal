#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep
from romancal.wfi18_transient.wfi18_transient import correct_anomaly
from romancal.stpipe.utilities import record_step_status

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

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        if self.save_results:
            try:
                self.suffix = "wfi18_transient"
            except AttributeError:
                self["suffix"] = "wfi18_transient"

        # Ignore any data that is not WFI18
        if input_model.meta.instrument.detector != "WFI18":
            record_step_status(input_model, "wfi18_transient", status="N/A")
            return input_model

        log.info("Correcting the first read transient anomaly for WFI18")
        correct_anomaly(input_model, mask_rows=self.mask_rows)
        record_step_status(input_model, "wfi18_transient", status="COMPLETE")

        return input_model
