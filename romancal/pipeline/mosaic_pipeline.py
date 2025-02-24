#!/usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import romancal.datamodels.filetype as filetype
from romancal.datamodels import ModelLibrary

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
log = logging.getLogger()
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
    """

    # Define aliases to steps
    step_defs: ClassVar = {
        "flux": FluxStep,
        "skymatch": SkyMatchStep,
        "outlier_detection": OutlierDetectionStep,
        "resample": ResampleStep,
        "sourcecatalog": SourceCatalogStep,
    }

    # start the actual processing
    def process(self, input):
        """Process the Roman WFI data from Level 2 to Level 3"""

        log.info("Starting Roman mosaic level calibration pipeline ...")

        # open the input file
        file_type = filetype.check(input)
        if file_type == "asdf":
            raise TypeError("The level three pipeline input needs to be an association")

        if file_type == "asn":
            input = ModelLibrary(input, on_disk=self.on_disk)
            self.flux.suffix = "flux"
            result = self.flux.run(input)
            self.skymatch.suffix = "skymatch"
            result = self.skymatch.run(result)
            self.outlier_detection.suffix = "outlier_detection"
            result = self.outlier_detection.run(result)
            self.resample.suffix = "coadd"
            self.output_file = input.asn["products"][0]["name"]
            result = self.resample.run(result)
            self.sourcecatalog.output_file = self.output_file
            self.sourcecatalog.run(result)
            self.suffix = "coadd"
        # FIXME fails for other file_type results
        return result
