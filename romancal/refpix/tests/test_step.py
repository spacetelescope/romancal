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
    assert (result.data == regression.data).all()
    # Check the amp33
    assert (result.amp33 == regression.amp33).all()
    # Check left ref pix
    assert (
        result.border_ref_pix_left == regression.border_ref_pix_left
    ).all()
    # Check right ref pix
    assert (
        result.border_ref_pix_right == regression.border_ref_pix_right
    ).all()
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
