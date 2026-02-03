#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep
from romancal.dark_decay.dark_decay import apply_dark_decay

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["DarkDecayStep"]

log = logging.getLogger(__name__)


class DarkDecayStep(RomanStep):
    """
    Apply dark decay correction.
    """

    class_alias = "dark_decay"
    reference_file_types: ClassVar = []

    spec = """
    """

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        if self.save_results:
            try:
                self.suffix = "dark_decay"
            except AttributeError:
                self["suffix"] = "dark_decay"

        apply_dark_decay(input_model)
        input_model.meta.cal_step.dark_decay = "COMPLETE"

        return input_model
