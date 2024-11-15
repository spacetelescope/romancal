#!/usr/bin/env python
import logging

import romancal.datamodels.filetype as filetype

# step imports
from romancal.assign_wcs import AssignWcsStep
from romancal.dark_current import DarkCurrentStep
from romancal.datamodels.library import ModelLibrary
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
    and source_catalog.
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
        # make sure we update source catalog coordinates afer running TweakRegStep
        self.tweakreg.update_source_catalog_coordinates = True
        # make output filenames based on input filenames
        self.output_use_model = True

        log.info("Starting Roman exposure calibration pipeline ...")

        # determine the input type
        file_type = filetype.check(input)
        if file_type == "ModelLibrary":
            lib = input
        elif file_type == "asn":
            lib = ModelLibrary(input)
        else:
            lib = ModelLibrary([input])

        with lib:
            for model_index, model in enumerate(lib):
                self.dq_init.suffix = "dq_init"
                result = self.dq_init.run(model)

                del model

                result = self.saturation.run(result)

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
                        "source_catalog",
                        "dark",
                        "refpix",
                        "linearity",
                        "ramp_fit",
                        "jump",
                        "tweakreg",
                    ]:
                        result.meta.cal_step[step_str] = "SKIPPED"
                else:
                    result = self.refpix.run(result)
                    result = self.linearity.run(result)
                    result = self.dark_current.run(result)
                    result = self.rampfit.run(result)
                    result = self.assign_wcs.run(result)

                    if result.meta.exposure.type == "WFI_IMAGE":
                        result = self.flatfield.run(result)
                        result = self.photom.run(result)
                        result = self.source_catalog.run(result)
                    else:
                        log.info("Flat Field step is being SKIPPED")
                        log.info("Photom step is being SKIPPED")
                        log.info("Source Detection step is being SKIPPED")
                        log.info("Tweakreg step is being SKIPPED")
                        result.meta.cal_step.flat_field = "SKIPPED"
                        result.meta.cal_step.photom = "SKIPPED"
                        result.meta.cal_step.source_catalog = "SKIPPED"
                        result.meta.cal_step.tweakreg = "SKIPPED"

                lib.shelve(result, model_index)

        # Now that all the exposures are collated, run tweakreg
        # Note: this does not cover the case where the asn mixes imaging and spectral
        #          observations. This should not occur on-prem
        self.tweakreg.run(lib)

        log.info("Roman exposure calibration pipeline ending...")

        return lib
