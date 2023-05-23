from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from romancal.stpipe import RomanStep
    from roman_datamodels.datamodels import RampModel, RefpixRefModel

from .data import Coefficients, StandardView

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


@dataclass
class Control:
    """
    Convenience class for reading user configuration from the Step itself.

    - remove_offset: Turn on removing the data offset prior to the correction
    - remove_trends: Remove the time dependent boundary trends
    - cosine_interpolate: Interpolate the reference pixels using a cosine
                          interpolation.
    - fft_interpolate: Interpolate the padding on the reference pixels using
                       an apodized FFT interpolation.
    """

    remove_offset: bool = True
    remove_trends: bool = True
    cosine_interpolate: bool = True
    fft_interpolate: bool = False

    @classmethod
    def from_step(cls, step: RomanStep) -> Control:
        return cls(
            step.remove_offset,
            step.remove_trends,
            step.cosine_interpolate,
            step.fft_interpolate,
        )


def run_steps(
    datamodel: RampModel, refs: RefpixRefModel, control: Control
) -> RampModel:
    """
    Organize the steps to run the reference pixel correction.
    """

    # Read in the data from the datamodels
    log.info("Reading data from datamodel into single array")
    coeffs = Coefficients.from_ref(refs)
    standard = StandardView.from_datamodel(datamodel)

    # Remove offset from the data
    if control.remove_offset:
        standard = standard.remove_offset()
        log.info("Removed the general offset from data, to be re-applied later.")

    # Convert to channel view
    channel = standard.channels

    # Remove the boundary trends
    if control.remove_trends:
        channel = channel.remove_trends()
        log.info("Removed boundary trends (in time) from data.")

    # Cosine interpolate the the data
    if control.cosine_interpolate:
        channel = channel.cosine_interpolate()
        log.info("Cosine interpolated the reference pixels.")

    # FFT interpolate the data
    if control.fft_interpolate:
        channel = channel.fft_interpolate()
        log.info("FFT interpolated the reference pixel pads.")

    # Perform the reference pixel correction
    standard = channel.apply_correction(coeffs)
    log.info("Applied reference pixel correction")

    # Re-apply the offset (if necessary)
    standard.apply_offset()
    log.info("Re-applied the general offset (if removed) to the data.")

    # Write the data back to the datamodel
    standard.update(datamodel)
    log.info("Updated the datamodel with the corrected data.")

    return datamodel
