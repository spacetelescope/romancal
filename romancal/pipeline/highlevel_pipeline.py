#!/usr/bin/env python
import logging
from os.path import basename

import numpy as np
from roman_datamodels import datamodels as rdm

import romancal.datamodels.filetype as filetype

# step imports
from romancal.skymatch import SkyMatchStep
from romancal.outlier_detection import OutlierDetectionStep
from romancal.resample import ResampleStep
from romancal.lib import dqflags
from romancal.datamodels import ModelContainer


from ..stpipe import RomanPipeline

__all__ = ["HighLevelPipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class HighLevelPipeline(RomanPipeline):
    """
    HighLevelPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    skymatch, Outlierdetectionn and resample.
    """

    class_alias = "roman_hlp"

    spec = """
        save_results = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {
        "skymatch": SkyMatchStep,
        "outlierdet": OutlierDetectionStep,
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
        asn = None
        if file_type == "asdf":
            log.info("The level three pipeline input needs to be an association")
            return

        if file_type == "asn":
            asn = ModelContainer.read_asn(input)
            self.skymatch.suffix = "skymatch"
            result = self.skymatch(input)
            self.skymatch.suffix = "outlierdetection"
            result = self.outlierdetection(asn)
            self.skymatch.suffix = "i2d"
            result = self.resample(result)
            if input_filename:
                result.meta.filename = input_filename            
   
        return result

    def setup_output(self, input):
        """Determine the proper file name suffix to use later"""
        if input.meta.cal_step.ramp_fit == "COMPLETE":
            self.suffix = "cal"
            input.meta.filename = input.meta.filename.replace("uncal", self.suffix)
            input["output_file"] = input.meta.filename
            self.output_file = input.meta.filename
        else:
            self.suffix = "cal"

