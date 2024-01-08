#!/usr/bin/env python
import logging
from os.path import basename

import romancal.datamodels.filetype as filetype
from romancal.outlier_detection import OutlierDetectionStep
from romancal.resample import ResampleStep

# step imports
from romancal.skymatch import SkyMatchStep

from ..stpipe import RomanPipeline

__all__ = ["HighLevelPipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class HighLevelPipeline(RomanPipeline):
    """
    HighLevelPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    ``skymatch``, ``outlier_detection`` and ``resample``.
    """

    class_alias = "roman_hlp"

    spec = """
        save_results = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {
        "skymatch": SkyMatchStep,
        "outlier_detection": OutlierDetectionStep,
        "resample": ResampleStep,
    }

    # start the actual processing
    def process(self, input):
        """Process the Roman WFI data from Level 2 to Level 3"""

        log.info("Starting Roman high level calibration pipeline ...")
        if isinstance(input, str):
            input_filename = basename(input)
        else:
            input_filename = None

        # open the input file
        file_type = filetype.check(input)
        if file_type == "asdf":
            log.info("The level three pipeline input needs to be an association")
            return

        if file_type == "asn":
            self.skymatch.suffix = "skymatch"
            result = self.skymatch(input)
            self.skymatch.suffix = "outlier_detection"
            # result = self.outlier_detection(input)
            self.skymatch.suffix = "i2d"
            result = self.resample(input)
            if input_filename:
                result.meta.filename = input_filename

        return result
