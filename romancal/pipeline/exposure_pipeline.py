#!/usr/bin/env python
import logging
from os.path import basename

import numpy as np
from roman_datamodels import datamodels as rdm
from roman_datamodels.dqflags import group

import romancal.datamodels.filetype as filetype

# step imports
from romancal.assign_wcs import AssignWcsStep
from romancal.associations.asn_from_list import asn_from_list
from romancal.associations.load_asn import load_asn
from romancal.dark_current import DarkCurrentStep
from romancal.dq_init import dq_init_step
from romancal.flatfield import FlatFieldStep
from romancal.lib.basic_utils import is_fully_saturated
from romancal.linearity import LinearityStep
from romancal.photom import PhotomStep
from romancal.ramp_fitting import ramp_fit_step
from romancal.refpix import RefPixStep
from romancal.saturation import SaturationStep
from romancal.source_catalog import SourceCatalogStep
from romancal.tweakreg import TweakRegStep

from ..stpipe import RomanPipeline

__all__ = ["ExposurePipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class ExposurePipeline(RomanPipeline):
    """
    ExposurePipeline: Apply all calibration steps to raw Roman WFI
    ramps to produce a 2-D slope product. Included steps are:
    dq_init, saturation, linearity, dark current, jump detection, ramp_fit,
    assign_wcs, flatfield (only applied to WFI imaging data), photom,
    and source_detection.
    """

    class_alias = "roman_elp"

    spec = """
        save_calibrated_ramp = boolean(default=False)
        save_results = boolean(default=False)
        suffix = string(default="cal")
    """

    # Define aliases to steps
    step_defs = {
        "dq_init": dq_init_step.DQInitStep,
        "saturation": SaturationStep,
        "refpix": RefPixStep,
        "linearity": LinearityStep,
        "dark_current": DarkCurrentStep,
        "rampfit": ramp_fit_step.RampFitStep,
        "assign_wcs": AssignWcsStep,
        "flatfield": FlatFieldStep,
        "photom": PhotomStep,
        "source_catalog": SourceCatalogStep,
        "tweakreg": TweakRegStep,
    }

    # start the actual processing
    def process(self, input):
        """Process the Roman WFI data"""

        # make sure source_catalog returns the updated datamodel
        self.source_catalog.return_updated_model = True

        log.info("Starting Roman exposure calibration pipeline ...")
        if isinstance(input, str):
            input_filename = basename(input)
        else:
            input_filename = None

        # determine the input type
        file_type = filetype.check(input)
        if file_type == "asn":
            with open(input_filename) as f:
                asn = load_asn(f)
        elif file_type == "asdf":
            try:
                # set the product name based on the input filename
                asn = asn_from_list([input], product_name=input_filename.split(".")[0])
                file_type = "asn"
            except TypeError:
                log.debug("Error opening file:")
                return

        # Build a list of observations to process
        expos_file = []
        n_members = 0
        # extract the members from the asn to run the files through the steps
        results = []
        tweakreg_input = []

        # populate expos_file
        if file_type == "asn":
            for product in asn["products"]:
                # extract the members from the asn to run the files through the steps
                n_members = len(product["members"])
                for member in product["members"]:
                    expos_file.append(member["expname"])
        elif file_type == "DataModel":
            expos_file.append(input)

        for in_file in expos_file:
            if isinstance(in_file, str):
                input_filename = basename(in_file)
                # Open the file
                input = rdm.open(in_file)
            elif isinstance(in_file, rdm.DataModel):
                input_filename = in_file.meta.filename
            else:
                input_filename = None

            log.info(f"Processing a WFI exposure {input_filename}")

            self.dq_init.suffix = "dq_init"
            result = self.dq_init(input)

            if input_filename:
                result.meta.filename = input_filename

            result = self.saturation(result)

            # Test for fully saturated data
            if is_fully_saturated(result):
                # Return fully saturated image file (stopping pipeline)
                log.info("All pixels are saturated. Returning a zeroed-out image.")

                #    if is_fully_saturated(result):
                # Set all subsequent steps to skipped
                for step_str in [
                    "assign_wcs",
                    "flat_field",
                    "photom",
                    "source_detection",
                    "dark",
                    "refpix",
                    "linearity",
                    "ramp_fit",
                    "jump",
                    "tweakreg",
                ]:
                    result.meta.cal_step[step_str] = "SKIPPED"

                results.append(result)
                return results

            result = self.refpix(result)
            result = self.linearity(result)
            result = self.dark_current(result)
            result = self.rampfit(result)
            result = self.assign_wcs(result)

            if result.meta.exposure.type == "WFI_IMAGE":
                result = self.flatfield(result)
                result = self.photom(result)
                result = self.source_catalog(result)
                tweakreg_input.append(result)
                log.info(
                    f"Number of models to tweakreg:   {len(tweakreg_input), n_members}"
                )
            else:
                log.info("Flat Field step is being SKIPPED")
                log.info("Photom step is being SKIPPED")
                log.info("Source Detection step is being SKIPPED")
                log.info("Tweakreg step is being SKIPPED")
                result.meta.cal_step.flat_field = "SKIPPED"
                result.meta.cal_step.photom = "SKIPPED"
                result.meta.cal_step.source_detection = "SKIPPED"
                result.meta.cal_step.tweakreg = "SKIPPED"

            self.output_use_model = True
            results.append(result)

        # Now that all the exposures are collated, run tweakreg
        # Note: this does not cover the case where the asn mixes imaging and spectral
        #          observations. This should not occur on-prem
        result = self.tweakreg(results)

        log.info("Roman exposure calibration pipeline ending...")

        return results

    def create_fully_saturated_zeroed_image(self, input_model):
        """
        Create zeroed-out image file
        """
        # The set order is: data, dq, var_poisson, var_rnoise, err
        fully_saturated_model = ramp_fit_step.create_image_model(
            input_model,
            (
                np.zeros(input_model.data.shape[1:], dtype=input_model.data.dtype),
                input_model.pixeldq | input_model.groupdq[0] | group.SATURATED,
                np.zeros(input_model.err.shape[1:], dtype=input_model.err.dtype),
                np.zeros(input_model.err.shape[1:], dtype=input_model.err.dtype),
                np.zeros(input_model.err.shape[1:], dtype=input_model.err.dtype),
            ),
        )

        # Set all subsequent steps to skipped
        for step_str in [
            "linearity",
            "dark",
            "ramp_fit",
            "assign_wcs",
            "flat_field",
            "photom",
            "source_detection",
            "tweakreg",
        ]:
            fully_saturated_model.meta.cal_step[step_str] = "SKIPPED"

        # Return zeroed-out image file
        return fully_saturated_model
