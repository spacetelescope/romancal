import pytest

from romancal.refpix import RefPixStep
from romancal.refpix.refpix import run_steps


@pytest.mark.parametrize(
    "remove_offset, remove_trends, cosine_interpolate, fft_interpolate",
    [(True, True, True, False), (True, True, True, True)],
)
def test_refpix_step(
    datamodel,
    ref_pix_ref,
    remove_offset,
    remove_trends,
    cosine_interpolate,
    fft_interpolate,
):
    # Setup and run the correction outside of the step, this is
    # already regression tested.
    regression_datamodel = datamodel.copy()
    regression_ref_pix_ref = ref_pix_ref.copy()
    regression = run_steps(
        regression_datamodel,
        regression_ref_pix_ref,
        remove_offset,
        remove_trends,
        cosine_interpolate,
        fft_interpolate,
    )

    # Run the step
    result = RefPixStep.call(
        datamodel,
        override_refpix=ref_pix_ref,
        remove_offset=remove_offset,
        remove_trends=remove_trends,
        cosine_interpolate=cosine_interpolate,
        fft_interpolate=fft_interpolate,
    )

    # Check the outputs are distinct
    assert result is not regression

    # Check the data
    assert (result.data.value == regression.data.value).all()
    assert result.data.unit == regression.data.unit
    # Check the amp33
    assert (result.amp33.value == regression.amp33.value).all()
    assert result.amp33.unit == regression.amp33.unit
    # Check left ref pix
    assert (
        result.border_ref_pix_left.value == regression.border_ref_pix_left.value
    ).all()
    assert result.border_ref_pix_left.unit == regression.border_ref_pix_left.unit
    # Check right ref pix
    assert (
        result.border_ref_pix_right.value == regression.border_ref_pix_right.value
    ).all()
    assert result.border_ref_pix_right.unit == regression.border_ref_pix_right.unit
    #
    # Run the step with reffile = N/A
    result = RefPixStep.call(
        datamodel,
        override_refpix="N/A",
        remove_offset=remove_offset,
        remove_trends=remove_trends,
        cosine_interpolate=cosine_interpolate,
        fft_interpolate=fft_interpolate,
    )
    assert result.meta.cal_step.refpix == "SKIPPED"
