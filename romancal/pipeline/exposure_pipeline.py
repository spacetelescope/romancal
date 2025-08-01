#!/usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels.dqflags import group

import romancal.datamodels.filetype as filetype

# step imports
from romancal.assign_wcs import AssignWcsStep
from romancal.dark_current import DarkCurrentStep
from romancal.datamodels.library import ModelLibrary
from romancal.dq_init import dq_init_step
from romancal.flatfield import FlatFieldStep
from romancal.lib.basic_utils import is_fully_saturated
from romancal.lib.save_wcs import save_wfiwcs
from romancal.linearity import LinearityStep
from romancal.photom import PhotomStep
from romancal.ramp_fitting import ramp_fit_step
from romancal.refpix import RefPixStep
from romancal.saturation import SaturationStep
from romancal.source_catalog import SourceCatalogStep
from romancal.tweakreg import TweakRegStep

from ..stpipe import RomanPipeline

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["ExposurePipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ExposurePipeline(RomanPipeline):
    """
    ExposurePipeline: Apply all calibration steps to raw Roman WFI
    ramps to produce a 2-D slope product. Included steps are documented
    in the ``step_defs``.
    """

    class_alias = "roman_elp"

    spec = """
        save_l1_wcs = boolean(default=False)
        save_results = boolean(default=False)
        suffix = string(default="cal")
    """

    # Define aliases to steps
    step_defs: ClassVar = {
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
        # tweakreg currently holds responsibiility for creating the L1 WCS files.
        self.tweakreg.save_l1_wcs = True
        # make output filenames based on input filenames
        self.output_use_model = True

        log.info("Starting Roman exposure calibration pipeline ...")

        # determine the input type
        file_type = filetype.check(input)
        return_lib = True
        if file_type == "ModelLibrary":
            lib = input
        elif file_type == "asn":
            lib = ModelLibrary(input)
        else:
            # for a non-asn non-library input process it as a library
            lib = ModelLibrary([input])
            # but return it as a datamodel
            return_lib = False

        # Flag to track if any of the input models are fully saturated
        any_saturated = False

        with lib:
            for model_index, model in enumerate(lib):
                self.dq_init.suffix = "dq_init"
                result = self.dq_init.run(model)

                del model

                result = self.saturation.run(result)

                if is_fully_saturated(result):
                    log.info("All pixels are saturated. Returning a zeroed-out image.")
                    result = self.create_fully_saturated_zeroed_image(result)

                    # Track that we've seen a fully saturated input
                    any_saturated = True
                    log.warning(
                        "tweakreg will not be run due to a fully saturated input"
                    )
                else:
                    result = self.refpix.run(result)
                    result = self.linearity.run(result)
                    result = self.rampfit.run(result)
                    result = self.dark_current.run(result)
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

                if any_saturated:
                    # the input association contains a fully saturated model
                    # where source_catalog can't be run which means we
                    # also can't run tweakreg.
                    result.meta.cal_step.tweakreg = "SKIPPED"
                lib.shelve(result, model_index)

        # Now that all the exposures are collated, run tweakreg
        # Note: this does not cover the case where the asn mixes imaging and spectral
        #          observations. This should not occur on-prem
        if not any_saturated:
            self.tweakreg.run(lib)

        # Write out the WfiWcs products
        if self.save_l1_wcs:
            save_wfiwcs(self, lib)

        log.info("Roman exposure calibration pipeline ending...")

        # return a ModelLibrary
        if return_lib:
            return lib

        # or a DataModel (for non-asn non-lib inputs)
        with lib:
            model = lib.borrow(0)
            lib.shelve(model)
        return model

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
            "refpix",
            "linearity",
            "dark",
            "ramp_fit",
            "assign_wcs",
            "flat_field",
            "photom",
            "source_catalog",
        ]:
            fully_saturated_model.meta.cal_step[step_str] = "SKIPPED"

        # Return zeroed-out image file
        return fully_saturated_model
