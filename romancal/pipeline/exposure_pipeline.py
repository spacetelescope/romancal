#!/usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import roman_datamodels.datamodels as rdm
from roman_datamodels.dqflags import group, pixel

# step imports
from romancal.assign_wcs import AssignWcsStep
from romancal.dark_current import DarkCurrentStep
from romancal.dark_decay import DarkDecayStep
from romancal.datamodels.fileio import open_dataset
from romancal.datamodels.library import ModelLibrary
from romancal.dq_init import dq_init_step
from romancal.flatfield import FlatFieldStep
from romancal.lib.save_wcs import save_wfiwcs
from romancal.linearity import LinearityStep
from romancal.photom import PhotomStep
from romancal.ramp_fitting import ramp_fit_step
from romancal.refpix import RefPixStep
from romancal.saturation import SaturationStep
from romancal.source_catalog import SourceCatalogStep
from romancal.tweakreg import TweakRegStep
from romancal.wfi18_transient import WFI18TransientStep

from ..stpipe import RomanPipeline

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["ExposurePipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _is_fully_saturated(model):
    """
    Check to see if all data pixels are flagged as saturated.
    """

    if np.all(np.bitwise_and(model.groupdq, group.SATURATED) == group.SATURATED):
        return True
    elif np.all(np.bitwise_and(model.pixeldq, pixel.SATURATED) == pixel.SATURATED):
        return True

    return False


class ExposurePipeline(RomanPipeline):
    """
    ExposurePipeline: Apply all calibration steps to raw Roman WFI
    ramps to produce a 2-D slope product. Included steps are documented
    in the ``step_defs``.
    """

    class_alias = "roman_elp"

    spec = """
        save_results = boolean(default=False)
        on_disk = boolean(default=False)
        suffix = string(default="cal")
    """

    # Define aliases to steps
    step_defs: ClassVar = {
        "dq_init": dq_init_step.DQInitStep,
        "saturation": SaturationStep,
        "refpix": RefPixStep,
        "dark_decay": DarkDecayStep,
        "wfi18_transient": WFI18TransientStep,
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
    def process(self, dataset):
        """Process the Roman WFI data"""

        # make output filenames based on input filenames
        self.output_use_model = True

        log.info("Starting Roman exposure calibration pipeline ...")

        # determine the input type
        lib, input_type = open_dataset(
            dataset,
            update_version=self.update_version,
            return_type=True,
            as_library=True,
            open_kwargs={"on_disk": self.on_disk},
        )
        return_lib = input_type in ("ModelLibrary", "asn")

        catalogs = []
        segmentations = []

        with lib:
            for model_index, model in enumerate(lib):
                result, run_source_catalog = self._process_model(model)

                # now handle source_catalog
                if not run_source_catalog or self.source_catalog.skip:
                    # WFI_WFSC doesn't get a source catalog (and therefore also no tweakreg)
                    result.meta.cal_step.source_catalog = "SKIPPED"
                    catalog, segmentation = None, None
                else:
                    # WFI_IMAGE and WFI_LOLO get source catalog
                    catalog, segmentation = self.source_catalog.run(result)

                if not self.tweakreg.skip and catalog is not None:
                    # attach the catalog to the model so tweakreg can see it
                    if "source_catalog" not in result.meta:
                        result.meta["source_catalog"] = {}
                    result.meta.source_catalog.tweakreg_catalog = catalog.source_catalog

                lib.shelve(result, model_index)
                catalogs.append(catalog)
                segmentations.append(segmentation)

        # Now that all the exposures are collated, run tweakreg
        self.tweakreg.run(lib)

        # tweakreg was run, update catalog positions
        with lib:
            for model_index, model in enumerate(lib):
                if model.meta.cal_step.tweakreg == "COMPLETE":
                    catalog = catalogs[model_index]
                    if catalog is not None:
                        self.tweakreg._update_catalog_coordinates(
                            catalog.source_catalog, model.meta.wcs
                        )
                        # record the name of the catalog if it is going to be saved
                        if self.save_results:
                            catalog_filename = self.make_output_path(
                                catalog.meta.filename, suffix="cat", ext="parquet"
                            )
                            model.meta.source_catalog.tweakreg_catalog_name = (
                                catalog_filename
                            )
                lib.shelve(model)

        log.info("Roman exposure calibration pipeline ending...")

        # return a ModelLibrary
        if return_lib:
            return (
                lib,
                ModelLibrary([c for c in catalogs if c is not None]),
                ModelLibrary([s for s in segmentations if s is not None]),
            )

        # or a DataModel (for non-asn non-lib inputs)
        with lib:
            model = lib.borrow(0)
            catalog = catalogs[0]
            segmentation = segmentations[0]
            lib.shelve(model, modify=False)
        return model, catalog, segmentation

    def _process_model(self, model):
        """
        Run all per-model calibration steps.

        Returns the model and a boolean indicating if source catalog should be run.
        """
        self.dq_init.suffix = "dq_init"
        result = self.dq_init.run(model)
        if model is not result:
            # dq_init converted this to a new model type so close the input
            model.close()
            del model

        result = self.saturation.run(result)

        if _is_fully_saturated(result):
            log.info("All pixels are saturated. Returning a zeroed-out image.")
            return self.create_fully_saturated_zeroed_image(result), False

        result = self.refpix.run(result)
        result = self.dark_decay.run(result)
        result = self.wfi18_transient.run(result)
        result = self.linearity.run(result)
        result = self.rampfit.run(result)
        result = self.dark_current.run(result)
        result = self.assign_wcs.run(result)
        result = self.photom.run(result)

        # WFI_FLAT, WFI_SPECTRAL, WFI_IM_DARK, WFI_SP_DARK stop here
        if result.meta.exposure.type not in ("WFI_IMAGE", "WFI_LOLO", "WFI_WFSC"):
            result.meta.cal_step.flat_field = "SKIPPED"
            result.meta.cal_step.source_catalog = "SKIPPED"
            return result, False

        return self.flatfield.run(result), result.meta.exposure.type in (
            "WFI_IMAGE",
            "WFI_LOLO",
        )

    def save_model(self, model, **kwargs):
        suffix = kwargs.get("suffix", None)
        # depending on model set suffix and ext
        if isinstance(
            model,
            (
                rdm.ForcedImageSourceCatalogModel,
                rdm.ImageSourceCatalogModel,
            ),
        ):
            kwargs["ext"] = "parquet"
            if suffix is None:
                suffix = "cat"
        elif isinstance(model, rdm.SegmentationMapModel):
            if suffix is None:
                suffix = "segm"
            kwargs["suffix"] = kwargs.get("suffix", "segm")
        elif isinstance(model, rdm.ImageModel):
            save_wfiwcs(self, model, force=True)
            if suffix is None:
                suffix = self.suffix

        kwargs["suffix"] = suffix

        # strip the index since these all have different extensions
        kwargs.pop("idx", None)

        return super().save_model(model, **kwargs)

    def create_fully_saturated_zeroed_image(self, input_model):
        """
        Create zeroed-out image file
        """
        # Make a throw-away model to get the expected datatypes
        fake_model = rdm.ImageModel.create_fake_data()

        # Create a dictionary for fully saturated data
        slopes = np.zeros(input_model.data.shape[1:], dtype=fake_model.data.dtype)
        dq = input_model.pixeldq | input_model.groupdq[0] | group.SATURATED
        err = np.zeros(input_model.data.shape[1:], dtype=fake_model.err.dtype)
        image_info_allsat = {
            "slope": slopes,
            "dq": dq,
            "var_poisson": err,
            "var_rnoise": err,
            "err": err,
        }

        fully_saturated_model = ramp_fit_step._create_image_model(
            input_model, image_info_allsat
        )

        # Set all subsequent steps to skipped
        for step_str in [
            "refpix",
            "dark_decay",
            "wfi18_transient",
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
