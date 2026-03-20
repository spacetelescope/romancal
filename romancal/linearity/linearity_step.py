"""
Apply linearity correction to a science image
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from roman_datamodels import datamodels as rdd
from roman_datamodels.dqflags import group, pixel
from stcal.linearity.linearity import (
    apply_polynomial,
    linearity_correction,
    prepare_coefficients,
)

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


def _linearity_correction_lowmem(
    data,
    gdq,
    pdq,
    lin_coeffs,
    lin_dq,
    dqflags,
    ilin_coeffs,
    additional_correction=None,
    read_pattern=None,
):
    """Memory-efficient read-level linearity correction.

    Drop-in replacement for stcal's ``linearity_correction`` that processes
    reads one at a time instead of materializing all reads per resultant
    simultaneously.  This reduces peak temporaries from
    ~3 * n_reads * nrows * ncols * 4 bytes to ~3 * nrows * ncols * 4 bytes.

    Two-pass algorithm per resultant:
      Pass 1 – accumulate mean(ilin_poly(predicted_read)) to compute offset
      Pass 2 – for each read, apply offset → optional INL → optional sat cap
               → lin_poly → accumulate corrected sum

    Parameters match ``stcal.linearity.linearity.linearity_correction``
    (subset used by Roman: single integration, no zeroframe).
    """
    # Prepare coefficients (handles NaN, zero, NO_LIN_CORR flags)
    lin_coeffs, new_pdq = prepare_coefficients(lin_coeffs, lin_dq, pdq, dqflags)
    ilin_coeffs, new_pdq = prepare_coefficients(ilin_coeffs, lin_dq, new_pdq, dqflags)

    nints = data.shape[0]
    for integration in range(nints):
        _linearity_int_lowmem(
            data[integration],
            gdq[integration],
            lin_coeffs,
            dqflags,
            ilin_coeffs=ilin_coeffs,
            additional_correction=additional_correction,
            read_pattern=read_pattern,
        )

    return data, new_pdq, None


def _linearity_int_lowmem(
    data,
    gdq,
    lin_coeffs,
    dqflags,
    ilin_coeffs,
    additional_correction=None,
    read_pattern=None,
):
    """Process one integration with incremental read-level correction.

    Parameters
    ----------
    data : ndarray, shape (ngroups, nrows, ncols)
    gdq : ndarray, shape (ngroups, nrows, ncols)
    lin_coeffs, ilin_coeffs : ndarray, shape (ncoeffs, nrows, ncols)
    dqflags : dict
    additional_correction : callable or None
    read_pattern : list of lists
    """
    ngroups, nrows, ncols = data.shape
    mean_read_resultant = np.array([np.mean(reads) for reads in read_pattern])

    # With fewer than 2 resultants we cannot estimate a count rate,
    # so fall back to simple polynomial linearization.
    if ngroups < 2:
        for i in range(ngroups):
            corrected = apply_polynomial(data[i], lin_coeffs, gdq[i], dqflags)
            resultant_saturated = (gdq[i] & dqflags["SATURATED"]) != 0
            data[i] = np.where(resultant_saturated, data[i], corrected)
        return data

    # Identify saturated pixels and count unsaturated resultants
    is_saturated = (gdq & dqflags["SATURATED"]) != 0
    n_unsaturated = np.sum(~is_saturated, axis=0) - 1
    n_unsaturated = np.clip(n_unsaturated, 2, ngroups)

    # Linearize first resultant and compute count rate
    firstread_lin = apply_polynomial(data[0], lin_coeffs, gdq[0], dqflags)

    lastvalidresultant = n_unsaturated - 1
    iy, ix = np.meshgrid(np.arange(nrows), np.arange(ncols), indexing="ij")
    last_idx = (lastvalidresultant[iy, ix], iy, ix)

    lastread_lin = apply_polynomial(data[last_idx], lin_coeffs, gdq[last_idx], dqflags)
    d_reads = mean_read_resultant[last_idx[0]] - mean_read_resultant[0]
    countrate = (lastread_lin - firstread_lin) / d_reads

    del lastread_lin, d_reads, last_idx, lastvalidresultant, n_unsaturated, is_saturated

    # Scratch arrays reused across resultants
    read_buf = np.empty((nrows, ncols), dtype=data.dtype)

    for i in range(ngroups):
        reads_in_group = read_pattern[i]
        n_reads = len(reads_in_group)
        reads_since_first = np.array(reads_in_group) - mean_read_resultant[0]

        # --- Pass 1: compute mean of unlinearized reads to get offset ---
        unlin_sum = np.zeros((nrows, ncols), dtype=np.float64)
        for dt in reads_since_first:
            np.add(firstread_lin, countrate * dt, out=read_buf)
            unlin_sum += apply_polynomial(read_buf[np.newaxis], ilin_coeffs)[0]
        predicted_cts = unlin_sum / n_reads
        offset = data[i] - predicted_cts
        del unlin_sum, predicted_cts

        # --- Pass 2: apply offset, corrections, linearize, accumulate ---
        corrected_sum = np.zeros((nrows, ncols), dtype=np.float64)
        for dt in reads_since_first:
            np.add(firstread_lin, countrate * dt, out=read_buf)
            read_unlin = apply_polynomial(read_buf[np.newaxis], ilin_coeffs)[0]
            read_unlin += offset

            if additional_correction is not None:
                read_unlin += additional_correction(read_unlin[np.newaxis])[0]

            corrected_sum += apply_polynomial(read_unlin, lin_coeffs)
            del read_unlin

        corrected_resultant = (corrected_sum / n_reads).astype(data.dtype)
        del corrected_sum

        resultant_saturated = (gdq[i] & dqflags["SATURATED"]) != 0
        data[i] = np.where(resultant_saturated, data[i], corrected_resultant)

    return data


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

            # Use memory-efficient implementation when inverse linearity
            # coefficients are available (read-level correction).  The stcal
            # version materializes all reads per resultant simultaneously,
            # which can exceed 6 GB for 32-read resultants on a 4096×4096
            # detector.  Our implementation processes reads one at a time.
            if ilin_coeffs is not None:
                new_data, new_pdq, _ = _linearity_correction_lowmem(
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
            else:
                new_data, new_pdq, _ = linearity_correction(
                    input_model.data,
                    gdq,
                    pdq,
                    lin_coeffs,
                    lin_dq,
                    pixel,
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
        input_model.data[m] = np.clip(input_model.data[m], -1e6, 1e6)
        input_model.groupdq[m] |= group.DO_NOT_USE
        nbad = np.sum(m)
        log.warning(f"Flagged {nbad} spurious values outside remotely plausible range.")

        # Update the step status
        input_model.meta.cal_step["linearity"] = "COMPLETE"

        if self.save_results:
            try:
                self.suffix = "linearity"
            except AttributeError:
                self["suffix"] = "linearity"
        return input_model
