#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import roman_datamodels as rdm

from romancal.datamodels.fileio import open_dataset
from romancal.photom import photom
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["PhotomStep"]

log = logging.getLogger(__name__)


class PhotomStep(RomanStep):
    """
    PhotomStep: Module for loading photometric conversion information from
        reference files and attaching to the input science data model
    """

    class_alias = "photom"

    reference_file_types: ClassVar = ["photom"]

    def process(self, dataset):
        """Perform the photometric calibration step

        Parameters
        ----------
        dataset : Roman level 2 image datamodel (wfi_image-1.x.x)
            input roman datamodel

        Returns
        -------
        output_model : Roman level 2 image datamodel (wfi_image-1.x.x)
            output roman datamodel
        """

        input_model = open_dataset(dataset, update_version=self.update_version)

        # Get reference file
        reffile = self.get_reference_file(input_model, "photom")

        # Open the relevant photom reference file as a datamodel
        if reffile is not None and reffile != "N/A":
            # If there is a reference file, perform photom application
            photom_model = rdm.open(reffile)
            log.debug(f"Using PHOTOM ref file: {reffile}")

            # Do the correction
            if input_model.meta.exposure.type == "WFI_IMAGE":
                photom.apply_photom(input_model, photom_model)
                input_model.meta.cal_step.photom = "COMPLETE"
            else:
                log.warning("No photometric corrections for spectral data")
                log.warning("Photom step will be skipped")
                input_model.meta.cal_step.photom = "SKIPPED"
                input_model.meta.photometry.pixel_area = None
                input_model.meta.photometry.conversion_megajanskys = None
                input_model.meta.photometry.conversion_megajanskys_uncertainty = None
            photom_model.close()

        else:
            # Skip Photom step if no photom file
            log.warning("No PHOTOM reference file found")
            log.warning("Photom step will be skipped")
            input_model.meta.cal_step.photom = "SKIPPED"

        if self.save_results:
            try:
                self.suffix = "photom"
            except AttributeError:
                self["suffix"] = "photom"

        return input_model
