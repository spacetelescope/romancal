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
        # Because we're processing raw files, let's open without any
        # laziness; we need to propagate all of the bits into the ramps
        # anyway.  It also avoids bugs like:
        # https://github.com/spacetelescope/roman_datamodels/issues/631
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
                result = self._process_model(model)

                # now handle source_catalog
                if (
                    result.meta.exposure.type not in ("WFI_IMAGE", "WFI_LOLO")
                    or self.source_catalog.skip
                ):
                    # WFI_WFSC doesn't get a source catalog (and therefore also no tweakreg)
                    result.meta.cal_step.source_catalog = "SKIPPED"
                    catalog, segmentation = (
                        self.source_catalog._make_catalog_and_segmentation_models(
                            result
                        )
                    )
                    catalog.source_catalog = catalog.create_empty_catalog()
                    segmentation.data = np.zeros(model.data.shape[-2:], dtype=np.uint32)
                else:
                    # WFI_IMAGE and WFI_LOLO get source catalog
                    catalog, segmentation = self.source_catalog.run(result)

                if not self.tweakreg.skip:
                    # attach the catalog to the model so tweakreg can see it
                    if "source_catalog" not in model.meta:
                        result.meta["source_catalog"] = {}
                    result.meta.source_catalog.tweakreg_catalog = catalog.source_catalog

                lib.shelve(result, model_index)
                catalogs.append(catalog)
                segmentations.append(segmentation)

        # Now that all the exposures are collated, run tweakreg
        if not self.tweakreg.skip:
            self.tweakreg.run(lib)

            # tweakreg was run, update catalog positions
            with lib:
                for model_index, model in enumerate(lib):
                    catalog = catalogs[model_index]
                    self.tweakreg._update_catalog_coordinates(
                        catalog.source_catalog, model.meta.wcs
                    )
                    lib.shelve(model)

        log.info("Roman exposure calibration pipeline ending...")

        # return a ModelLibrary
        if return_lib:
            return lib, ModelLibrary(catalogs), ModelLibrary(segmentations)

        # or a DataModel (for non-asn non-lib inputs)
        with lib:
            model = lib.borrow(0)
            catalog = catalogs[0]
            segmentation = segmentations[0]
            lib.shelve(model, modify=False)
        return model, catalog, segmentation

    def _process_model(self, model):
        """Run all per-model calibration steps that return models and return the result."""
        self.dq_init.suffix = "dq_init"
        result = self.dq_init.run(model)
        if model is not result:
            # dq_init converted this to a new model type so close the input
            model.close()
            del model

        result = self.saturation.run(result)

        if _is_fully_saturated(result):
            log.info("All pixels are saturated. Returning a zeroed-out image.")
            return self.create_fully_saturated_zeroed_image(result)

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
            return result

        return self.flatfield.run(result)

    def save_model(self, model, **kwargs):
        # depending on model set suffix and ext
        if isinstance(
            model,
            (
                rdm.ForcedImageSourceCatalogModel,
                rdm.ImageSourceCatalogModel,
                rdm.ForcedMosaicSourceCatalogModel,
                rdm.MosaicSourceCatalogModel,
            ),
        ):
            kwargs["ext"] = "parquet"
            kwargs["suffix"] = kwargs.get("suffix", "cat")
        elif isinstance(
            model,
            (rdm.SegmentationMapModel, rdm.MosaicSegmentationMapModel),
        ):
            kwargs["suffix"] = kwargs.get("suffix", "segm")
        elif isinstance(model, rdm.ImageModel):
            save_wfiwcs(self, model, force=True)
            kwargs["suffix"] = kwargs.get("suffix", "cal")

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
