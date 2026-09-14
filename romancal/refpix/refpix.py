from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels.dqflags import pixel

if TYPE_CHECKING:
    from roman_datamodels.datamodels import RampModel, RefpixRefModel

from ._data import Coefficients, StandardView

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def run_steps(
    datamodel: RampModel,
    refs: RefpixRefModel,
    zero_bad_ref_pix: bool = True,
    remove_offset: bool = True,
    remove_trends: bool = True,
    cosine_interpolate: bool = True,
    fft_interpolate: bool = True,
) -> RampModel:
    """
    Organize the steps to run the reference pixel correction.

    Apart from the offset, the correction is independent between resultants, so
    the ramp is corrected one resultant at a time. The intermediate channel view
    and correction arrays are each larger than the ramp itself, so doing this
    keeps the peak memory usage proportional to a single resultant rather than
    to the whole ramp.
    """

    if zero_bad_ref_pix:
        _mask_bad_ref_pixels(datamodel)

    coeffs = Coefficients.from_ref(refs)

    # The offset is fit over all of the resultants, so it has to be accumulated
    # before any of them are corrected.
    offset = None
    if remove_offset:
        offset = StandardView.compute_offset(datamodel.data, datamodel.amp33)
        log.debug("Computed the general offset of the data, to be re-applied later.")

    log.debug(
        "Correcting the data one resultant at a time "
        f"(remove_trends={remove_trends}, cosine_interpolate={cosine_interpolate}, "
        f"fft_interpolate={fft_interpolate})."
    )

    for resultant in range(datamodel.data.shape[0]):
        # Read a single resultant from the datamodel into the standard view
        standard = StandardView.from_datamodel(datamodel, resultant)

        # Remove offset from the data
        if offset is not None:
            standard.remove_offset(offset)

        # Convert to channel view
        channel = standard.channels

        # Remove the boundary trends
        if remove_trends:
            channel = channel.remove_trends()

        # Cosine interpolate the the data
        if cosine_interpolate:
            channel = channel.cosine_interpolate()

        # FFT interpolate the data
        if fft_interpolate:
            channel = channel.fft_interpolate()

        # Perform the reference pixel correction
        standard = channel.apply_correction(coeffs)

        # Re-apply the offset (if necessary)
        standard.apply_offset()

        # Write the resultant back to the datamodel
        standard.update(datamodel, resultant)

    log.debug("Updated the datamodel with the corrected data.")

    return datamodel


def _mask_bad_ref_pixels(datamodel):
    """
    Zero out BAD_REF_PIXEL pixels before refpix correction.
    """
    dq = datamodel.pixeldq
    bad_ref_mask = (dq & pixel.BAD_REF_PIXEL) != 0

    # Protect against incorrectly shaped mock data in the unit tests.
    # Production datamodels enforce strict matching schemas.
    if bad_ref_mask.shape == datamodel.data.shape[-len(bad_ref_mask.shape) :]:
        datamodel.data[..., bad_ref_mask] = 0
    else:
        log.debug(
            "Skipping bad_ref_pixel masking: pixeldq shape %s does not match data trailing shape",
            bad_ref_mask.shape,
        )

    return datamodel
