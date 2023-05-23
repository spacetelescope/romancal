from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from romancal.stpipe import RomanStep
    from roman_datamodels.datamodels import RampModel, RefpixRefModel

from .data import Coefficients, StandardView


@dataclass
class Control:
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
    coeffs = Coefficients.from_ref(refs)
    standard = StandardView.from_datamodel(datamodel)

    # Remove offset from the data
    if control.remove_offset:
        standard = standard.remove_offset()

    # Convert to channel view
    channel = standard.channels

    # Remove the boundary trends
    if control.remove_trends:
        channel = channel.remove_trends()

    # Cosine interpolate the the data
    if control.cosine_interpolate:
        channel = channel.cosine_interpolate()

    # FFT interpolate the data
    if control.fft_interpolate:
        channel = channel.fft_interpolate()

    # Perform the reference pixel correction
    standard = channel.apply_correction(coeffs)

    # Re-apply the offset (if necessary)
    standard.apply_offset()

    # Write the data back to the datamodel
    standard.update(datamodel)

    return datamodel
