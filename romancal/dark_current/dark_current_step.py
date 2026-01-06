from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels import datamodels as rdm

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["DarkCurrentStep"]

log = logging.getLogger(__name__)


class DarkCurrentStep(RomanStep):
    """DarkCurrentStep: Performs dark current correction by subtracting
    dark current reference data from the input science data model.
    """

    class_alias = "dark"

    spec = """
        dark_output = output_file(default = None) # Dark corrected model
    """

    reference_file_types: ClassVar = ["dark"]

    def process(self, dataset):
        input_model = open_dataset(dataset)

        # Get the name of the dark reference file to use
        self.dark_name = self.get_reference_file(input_model, "dark")
        # Check for a valid reference file
        if self.dark_name == "N/A":
            log.warning("No DARK reference file found")
            log.warning("Dark current step will be skipped")
            input_model.meta.cal_step.dark = "SKIPPED"
            return input_model

        log.info("Using DARK reference file: %s", self.dark_name)

        # Open dark model
        with rdm.open(self.dark_name) as dark_model:
            # get the dark slope and dark slope error from the reference file & trim ref pixels
            dark_slope = dark_model.dark_slope[4:-4, 4:-4]
            dark_slope_err = dark_model.dark_slope_error[4:-4, 4:-4]

            # Do the dark correction
            input_model.data -= dark_slope
            input_model.err = np.sqrt(input_model.err**2 + (dark_slope_err) ** 2)
            input_model.dq |= dark_model.dq[4:-4, 4:-4]
            input_model.meta.cal_step.dark = "COMPLETE"

            # Save dark data to file
            if self.dark_output is not None:
                dark_model.save(self.dark_output)
                # not clear to me that this makes any sense for Roman

        if self.save_results:
            try:
                self.suffix = "darkcurrent"
            except AttributeError:
                self["suffix"] = "darkcurrent"

        return input_model
