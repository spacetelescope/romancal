"""
Apply linearity correction to a science image
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels import datamodels as rdd
from roman_datamodels.dqflags import pixel
from stcal.linearity.linearity import linearity_correction

from romancal.datamodels.fileio import open_dataset
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["LinearityStep"]

log = logging.getLogger(__name__)


def make_inl_correction(inl_model, ncols):
    """
    Create a callable for integral nonlinearity correction.

    Parameters
    ----------
    inl_model : datamodel
        The integral nonlinearity reference file model.
    ncols : int
        Number of columns in the data, used to determine which channels
        to extract.

    Returns
    -------
    callable
        A function that takes a 3D array (nreads, nrows, ncols) and returns
        a correction array of the same shape to be added to the data.
    """
    channel_width = 128
    lookup_values = inl_model.value.copy().astype("f4")
    channel_corrections = {}
    for start_col in range(0, ncols, channel_width):
        channel_num = start_col // channel_width + 1
        attr_name = f"science_channel_{channel_num:02d}"
        channel_corrections[channel_num] = getattr(
            inl_model.inl_table, attr_name
        ).correction.copy()

    def inl_correction(data):
        """Apply INL correction to data array."""
        result = np.zeros_like(data)
        for start_col in range(0, data.shape[-1], channel_width):
            channel_num = start_col // channel_width + 1
            correction = channel_corrections[channel_num]
            channel_data = data[..., start_col : start_col + channel_width]
            result[..., start_col : start_col + channel_width] = np.interp(
                channel_data, lookup_values, correction
            )
        return result

    return inl_correction


class LinearityStep(RomanStep):
    """
    LinearityStep: This step performs a correction for non-linear
    detector response, using the "classic" polynomial method.
    """

    class_alias = "linearity"

    reference_file_types: ClassVar = [
        "linearity",
        "inverselinearity",
        "integralnonlinearity",
    ]

    def process(self, dataset):
        input_model = open_dataset(dataset, update_version=self.update_version)

        # Get reference file names
        self.lin_name = self.get_reference_file(input_model, "linearity")
        self.ilin_name = self.get_reference_file(input_model, "inverselinearity")
        self.inl_name = self.get_reference_file(input_model, "integralnonlinearity")
        log.info("Using LINEARITY reference file: %s", self.lin_name)
        log.info("Using INVERSELINEARITY reference file: %s", self.ilin_name)
        log.info("Using INTEGRALNONLINEARITY reference file: %s", self.inl_name)

        # Check for valid reference files
        if self.lin_name == "N/A" or self.ilin_name == "N/A":
            log.warning("No LINEARITY or INVERSELINEARITY reference file found")
            log.warning("Linearity step will be skipped")
            input_model.meta.cal_step["linearity"] = "SKIPPED"
            return input_model

        # INL correction is optional
        inl_correction = None
        if self.inl_name != "N/A":
            with rdd.open(self.inl_name) as inl_model:
                inl_correction = make_inl_correction(
                    inl_model, input_model.data.shape[-1]
                )

        with (
            rdd.LinearityRefModel(self.lin_name, memmap=False) as lin_model,
            rdd.InverselinearityRefModel(self.ilin_name, memmap=False) as ilin_model,
        ):
            lin_coeffs = lin_model.coeffs
            lin_dq = lin_model.dq
            ilin_coeffs = ilin_model.coeffs
            read_pattern = input_model.meta.exposure.read_pattern

            gdq = input_model.groupdq[np.newaxis, :]
            pdq = input_model.pixeldq
            input_model.data = input_model.data[np.newaxis, :]

            # Call linearity correction function in stcal
            # The third return value is the processed zero frame which
            # Roman does not use.
            new_data, new_pdq, _ = linearity_correction(
                input_model.data,
                gdq,
                pdq,
                lin_coeffs,
                lin_dq,
                pixel,
                ilin_coeffs=ilin_coeffs,
                additional_correction=inl_correction,
                read_pattern=read_pattern,
            )

            input_model.data = new_data[0, :, :, :]
            input_model.pixeldq = new_pdq

        # FIXME: force all values in array to be at least vaguely sane.
        # This should not happen for good linearity corrections and linearity
        # correction flagging, but current reference files have issues that
        # cause more problems downstream.
        # Full well is 65k DN.  After linearity correction we can't be more than
        # a factor of several away from this.
        # Any points larger than 1e6 should be flagged.
        m = np.abs(input_model.data) > 1e6
        input_model.data[m] = np.sign(input_model.data[m]) * 1e6
        input_model.pixeldq[m] |= pixel.DO_NOT_USE

        # Update the step status
        input_model.meta.cal_step["linearity"] = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "linearity"
            except AttributeError:
                self["suffix"] = "linearity"
        return input_model
