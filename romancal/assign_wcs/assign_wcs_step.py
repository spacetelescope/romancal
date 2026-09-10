"""
Assign a gWCS object to a science image.

"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from romancal.datamodels.fileio import open_dataset

from ..stpipe import RomanStep
from .assign_wcs import load_wcs

if TYPE_CHECKING:
    from typing import ClassVar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["AssignWcsStep"]


class AssignWcsStep(RomanStep):
    """Assign a gWCS object to a science image."""

    class_alias = "assign_wcs"

    reference_file_types: ClassVar = ["distortion"]

    def process(self, dataset):
        reference_file_names = {}
        input_model = open_dataset(dataset, update_version=self.update_version)

        for reftype in self.reference_file_types:
            log.info(f"reftype, {reftype}")
            reffile = self.get_reference_file(input_model, reftype)
            # Check for a valid reference file
            if reffile == "N/A":
                log.warning("No DISTORTION reference file found")
                log.warning("Assign WCS step will be skipped")
                result = input_model.copy()
                result.meta.cal_step.assign_wcs = "SKIPPED"
                return result

            reference_file_names[reftype] = reffile if reffile else ""
        log.info("Using reference files: %s for assign_wcs", reference_file_names)
        result = load_wcs(input_model, reference_file_names)

        if self.save_results:
            try:
                self.suffix = "assignwcs"
            except AttributeError:
                self["suffix"] = "assignwcs"

        return result
