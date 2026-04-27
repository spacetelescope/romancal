#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import roman_datamodels as rdm

from romancal.datamodels.fileio import open_dataset
from romancal.dq_init.dq_initialization import do_dqinit
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["DQInitStep"]

log = logging.getLogger(__name__)


class DQInitStep(RomanStep):
    """Initialize the Data Quality extension from the
    mask reference file.

    The dq_init step initializes the pixeldq attribute of the input
    datamodel using the MASK reference file.  The dq attribute of the
    MASK model is bitwise OR'd with the pixeldq attribute of the input
    model.

    Also adjust data for data_encoding_offset and reference_read offsets,
    and mark guide-window affected pixels' masks.
    """

    class_alias = "dq_init"

    spec = """
        expand_gw_flagging = integer(default=0)  # expand guide window flagging this many pixels around guide window
    """

    reference_file_types: ClassVar = ["mask"]

    def process(self, dataset):
        """Perform the dq_init calibration step

        Parameters
        ----------
        dataset : Roman datamodel
            input roman datamodel

        Returns
        -------
        output_model : Roman datamodel
            result roman datamodel
        """
        # Open datamodel
        input_model = open_dataset(dataset, update_version=self.update_version)

        # Get reference file path
        reference_file_name = self.get_reference_file(input_model, "mask")
        if reference_file_name != "N/A" and reference_file_name is not None:
            # If there are mask files, perform dq step
            # Open the relevant reference files as datamodels
            reference_file_model = rdm.open(reference_file_name)
            log.debug(f"Using MASK ref file: {reference_file_name}")
        else:
            reference_file_model = None

        output_model = do_dqinit(
            input_model,
            mask=reference_file_model,
            expand_gw_flagging=self.expand_gw_flagging,
        )

        # Close the input and reference files
        input_model.close()

        try:
            reference_file_model.close()
        except AttributeError:
            pass

        if self.save_results:
            try:
                self.suffix = "dqinit"
            except AttributeError:
                self["suffix"] = "dqinit"

        return output_model
