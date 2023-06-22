import os
from dataclasses import asdict

import pytest

from romancal.refpix import RefpixStep
from romancal.refpix.refpix import Control, run_steps


@pytest.mark.bigdata
@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason=(
        "Roman CRDS servers are not currently available outside the internal network"
    ),
)
@pytest.mark.parametrize("control", [Control(), Control(fft_interpolate=True)])
def test_refpix_step(datamodel, ref_pix_ref, control):
    # Setup and run the correction outside of the step, this is
    # already regression tested.
    regression_datamodel = datamodel.copy()
    regression_ref_pix_ref = ref_pix_ref.copy()
    regression = run_steps(regression_datamodel, regression_ref_pix_ref, control)

    # Run the step
    result = RefpixStep.call(datamodel, override_refpix=ref_pix_ref, **asdict(control))

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
