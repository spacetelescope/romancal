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
    """

    # Read in the data from the datamodels
    log.debug("Reading data from datamodel into single array")

    if zero_bad_ref_pix:
        mask_bad_ref_pixels(datamodel)

    coeffs = Coefficients.from_ref(refs)
    standard = StandardView.from_datamodel(datamodel)

    # Remove offset from the data
    if remove_offset:
        standard = standard.remove_offset()
        log.debug("Removed the general offset from data, to be re-applied later.")

    # Convert to channel view
    channel = standard.channels

    # Remove the boundary trends
    if remove_trends:
        channel = channel.remove_trends()
        log.debug("Removed boundary trends (in time) from data.")

    # Cosine interpolate the the data
    if cosine_interpolate:
        channel = channel.cosine_interpolate()
        log.debug("Cosine interpolated the reference pixels.")

    # FFT interpolate the data
    if fft_interpolate:
        channel = channel.fft_interpolate()
        log.debug("FFT interpolated the reference pixel pads.")

    # Perform the reference pixel correction
    standard = channel.apply_correction(coeffs)
    log.debug("Applied reference pixel correction")

    # Re-apply the offset (if necessary)
    standard.apply_offset()
    log.debug("Re-applied the general offset (if removed) to the data.")

    # Write the data back to the datamodel
    standard.update(datamodel)
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
