#! /usr/bin/env python
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels import datamodels as rdd

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep
from romancal.dark_decay.dark_decay import subtract_dark_decay

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["DarkDecayStep"]

log = logging.getLogger(__name__)


class DarkDecayStep(RomanStep):
    """
    Apply dark decay correction.
    """

    class_alias = "dark_decay"
    reference_file_types: ClassVar = ["darkdecaysignal"]

    spec = """
    """

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        if self.save_results:
            try:
                self.suffix = "dark_decay"
            except AttributeError:
                self["suffix"] = "dark_decay"

        # Get the reference file
        reffile = self.get_reference_file(input_model, "darkdecaysignal")
        log.info("Using DARKDECAY reference file: %s", reffile)

        # Get the detector-specific decay table
        detector = input_model.meta.instrument.detector
        sca = int(detector[3:])
        with rdd.open(reffile) as reference:
            decay_table = getattr(reference.decay_table, detector)

        # Get exposure metadata
        frame_time = input_model.meta.exposure.frame_time
        read_pattern = input_model.meta.exposure.read_pattern

        subtract_dark_decay(
            input_model.data,
            decay_table.amplitude,
            decay_table.time_constant,
            frame_time,
            read_pattern,
            sca,
        )
        input_model.meta.cal_step.dark_decay = "COMPLETE"

        return input_model
