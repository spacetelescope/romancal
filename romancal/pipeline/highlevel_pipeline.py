#!/usr/bin/env python
import logging
from os.path import basename

import romancal.datamodels.filetype as filetype
from romancal.datamodels import ModelContainer

# step imports
from romancal.flux import FluxStep
from romancal.skymatch import SkyMatchStep
from romancal.outlier_detection import OutlierDetectionStep
from romancal.resample import ResampleStep

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
        "flux": FluxStep,
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
            input = ModelContainer(input)
            self.flux.suffix = "flux"
            result = self.flux(input)
            self.skymatch.suffix = "skymatch"
            result = self.skymatch(result)
            self.outlier_detection.suffix = "outlier_detection"
            result = self.outlier_detection(result)
            self.resample.suffix = "i2d"
            result = self.resample(result)
            self.suffix = "i2d"
            if input_filename:
                result.meta.filename = self.output_file
            self.output_file = input.asn_table['products'][0]['name']

        return result
