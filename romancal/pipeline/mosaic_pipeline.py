#!/usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels import datamodels as rdm

from romancal.datamodels.fileio import open_dataset

# step imports
from romancal.flux import FluxStep
from romancal.outlier_detection import OutlierDetectionStep
from romancal.resample import ResampleStep
from romancal.skymatch import SkyMatchStep
from romancal.source_catalog import SourceCatalogStep

from ..stpipe import RomanPipeline

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["MosaicPipeline"]

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class MosaicPipeline(RomanPipeline):
    """
    MosaicPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    ``flux``, ``skymatch``, ``outlier_detection``, ``resample`` and ``source catalog``.
    """

    class_alias = "roman_mos"
    spec = """
        save_results = boolean(default=False)
        on_disk = boolean(default=False)
        resample_on_skycell = boolean(default=True)
        suffix = string(default="coadd")
    """

    # Define aliases to steps
    step_defs: ClassVar = {
        "flux": FluxStep,
        "skymatch": SkyMatchStep,
        "outlier_detection": OutlierDetectionStep,
        "resample": ResampleStep,
        "source_catalog": SourceCatalogStep,
    }

    def save_model(self, model, **kwargs):
        if isinstance(
            model,
            (
                rdm.ForcedMosaicSourceCatalogModel,
                rdm.MosaicSourceCatalogModel,
            ),
        ):
            kwargs["ext"] = "parquet"
            kwargs["suffix"] = kwargs.get("suffix", "cat")
        elif isinstance(model, rdm.MosaicSegmentationMapModel):
            kwargs["suffix"] = kwargs.get("suffix", "segm")
        elif isinstance(model, rdm.MosaicModel):
            kwargs["suffix"] = kwargs.get("suffix", self.suffix)

        # strip the index since these all have different extensions
        kwargs.pop("idx", None)

        return super().save_model(model, **kwargs)

    # start the actual processing
    def process(self, dataset):
        """Process the Roman WFI data from Level 2 to Level 3"""

        log.info("Starting Roman mosaic level calibration pipeline ...")

        # open the input file
        library = open_dataset(
            dataset,
            update_version=self.update_version,
            as_library=True,
            open_kwargs={"on_disk": self.on_disk},
        )

        # propagate resample_on_skycell setting
        self.outlier_detection.resample_on_skycell = self.resample_on_skycell
        self.resample.resample_on_skycell = self.resample_on_skycell

        if self.save_results:
            # only set if save_results is True since setting
            # output_file will trigger results to be saved
            self.output_file = library.asn["products"][0]["name"]

        self.flux.suffix = "flux"
        result = self.flux.run(library)
        self.skymatch.suffix = "skymatch"
        result = self.skymatch.run(result)
        self.outlier_detection.suffix = "outlier_detection"
        result = self.outlier_detection.run(result)
        self.resample.suffix = "coadd"

        result = self.resample.run(result)
        catalog_and_segmentation = self.source_catalog.run(result)
        if self.source_catalog.skip:
            catalog, segmentation = None, None
        else:
            catalog, segmentation = catalog_and_segmentation
        return result, catalog, segmentation
