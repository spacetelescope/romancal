#!/usr/bin/env python
import logging
from os.path import basename

import numpy as np
from roman_datamodels import datamodels as rdm

import romancal.datamodels.filetype as filetype

# step imports
from romancal.skymatch import SkyMatchStep
from romancal.outlier_detection import OutlierDetectionStep
#from romancal.mosaic import MosaicStep
from romancal.lib import dqflags

from ..stpipe import RomanPipeline

__all__ = ["HighLevelPipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class HighLevelPipeline(RomanPipeline):
    """
    HighLevelPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    skymatch, Outlierdetectionn and mosaic.
    """

    class_alias = "roman_hlp"

    spec = """
        save_results = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {
        "skymatch": SkyMatchStep,
        "outlierdet": OutlierDetectionStep,
        #"moasic": MosaicStep,
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
            try:
                input = rdm.open(input)
            except TypeError:
                log.debug("Error opening file:")
                return

        if file_type == "asn":
            try:
                asn = LoadAsLevel2Asn.load(input, basename=self.output_file)
            except AssociationNotValidError:
                log.debug("Error opening file:")
                return

        # Build a list of observations to process
        expos_file = []
        if file_type == "asdf":
            expos_file = [input]
        elif file_type == "asn":
            for product in asn["products"]:
                for member in product["members"]:
                    expos_file.append(member["expname"])

        results = []
        for in_file in expos_file:
            if isinstance(in_file, str):
                input_filename = basename(in_file)
                log.info(f"Input file name: {input_filename}")
            else:
                input_filename = None

            # Open the file
            input = rdm.open(in_file)
            log.info(f"Processing a WFI exposure {in_file}")

            self.skymatch.suffix = "skymatch"
            result = self.skymatch(input)
            if input_filename:
                result.meta.filename = input_filename
            result = self.outlierdetection(result)

            result = self.mosaic(result)
            # setup output_file for saving
            self.setup_output(result)
            log.info("Roman exposure calibration pipeline ending...")

            self.output_use_model = True
            results.append(result)

        return results

    def setup_output(self, input):
        """Determine the proper file name suffix to use later"""
        if input.meta.cal_step.ramp_fit == "COMPLETE":
            self.suffix = "cal"
            input.meta.filename = input.meta.filename.replace("uncal", self.suffix)
            input["output_file"] = input.meta.filename
            self.output_file = input.meta.filename
        else:
            self.suffix = "ramp"

