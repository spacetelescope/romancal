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
    """

    # Define aliases to steps
    step_defs: ClassVar = {
        "flux": FluxStep,
        "skymatch": SkyMatchStep,
        "outlier_detection": OutlierDetectionStep,
        "resample": ResampleStep,
        "source_catalog": SourceCatalogStep,
    }

    # start the actual processing
    def process(self, input_data):
        """Process the Roman WFI data from Level 2 to Level 3"""

        log.info("Starting Roman mosaic level calibration pipeline ...")

        # open the input file
        file_type = filetype.check(input_data)
        if file_type == "asn":
            library = ModelLibrary(input_data, on_disk=self.on_disk)
        elif file_type == "ModelLibrary":
            library = input_data
        else:
            raise TypeError(
                "The level three pipeline input needs to be an association or ModelLibrary"
            )

        # propagate resample_on_skycell setting
        self.outlier_detection.resample_on_skycell = self.resample_on_skycell
        self.resample.resample_on_skycell = self.resample_on_skycell

        self.flux.suffix = "flux"
        result = self.flux.run(library)
        self.skymatch.suffix = "skymatch"
        result = self.skymatch.run(result)
        self.outlier_detection.suffix = "outlier_detection"
        result = self.outlier_detection.run(result)
        self.resample.suffix = "coadd"
        self.output_file = library.asn["products"][0]["name"]
        result = self.resample.run(result)
        self.source_catalog.output_file = self.output_file
        self.source_catalog.run(result)
        self.suffix = "coadd"
        return result
