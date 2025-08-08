"""
Flat-field a science image
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import roman_datamodels as rdm

from ..stpipe import RomanStep
from . import flat_field

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["FlatFieldStep"]

log = logging.getLogger(__name__)


class FlatFieldStep(RomanStep):
    """Flat-field a science image using a flatfield reference image."""

    class_alias = "flat_field"
    spec = """
        include_var_flat = boolean(default=False) # include flat field variance
    """

    reference_file_types: ClassVar = ["flat"]

    def process(self, input_model):
        if not isinstance(input_model, rdm.DataModel):
            input_model = rdm.open(input_model)

        reference_file_name = self.get_reference_file(input_model, "flat")

        # Check for a valid reference file
        if reference_file_name == "N/A":
            log.warning("No FLAT reference file found")
            log.warning("Flat Field step will be skipped")
            reference_file_name = None

        if reference_file_name is not None:
            reference_file_model = rdm.open(reference_file_name)
            log.debug(f"Using FLAT ref file: {reference_file_name}")
        else:
            reference_file_model = None
            log.debug("Using no FLAT ref file")

        # Do the flat-field correction
        output_model = flat_field.do_correction(
            input_model, reference_file_model, include_var_flat=self.include_var_flat
        )

        # Close reference file
        if reference_file_model is not None:
            reference_file_model.close()

        if self.save_results:
            try:
                self.suffix = "flat"
            except AttributeError:
                self["suffix"] = "flat"

        return output_model
