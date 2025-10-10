"""Ramp Fitting tests involving MultiAccum Tables"""

import numpy as np
import pytest

from romancal.ramp_fitting import RampFitStep

from .common import SIMPLE_RESULTANTS, make_data

SIMPLE_EXPECTED_DEFAULT = {
    "data": np.array(
        [[0.52631587, 0.52631587], [0.23026317, 0.7236843]], dtype=np.float32
    ),
    "err": np.array(
        [[0.24262409, 0.24262409], [0.16048454, 0.28450054]], dtype=np.float32
    ),
    "var_poisson": np.array(
        [[0.05886428, 0.05886428], [0.02575312, 0.08093839]], dtype=np.float32
    ),
    "var_rnoise": np.array(
        [[2.164128e-06, 2.164128e-06], [2.164128e-06, 2.164128e-06]], dtype=np.float32
    ),
}
SIMPLE_EXPECTED_GAIN = {
    "data": np.array([[0.526316, 0.526316], [0.230263, 0.701852]], dtype=np.float32),
    "err": np.array([[0.108513, 0.108513], [0.071783, 0.124624]], dtype=np.float32),
    "var_poisson": np.array(
        [[1.1772858e-02, 1.1772858e-02], [5.150624e-03, 1.55289e-02]], dtype=np.float32
    ),
    "var_rnoise": np.array(
        [[2.164127e-06, 2.164127e-06], [2.164127e-06, 2.190606e-06]], dtype=np.float32
    ),
}
SIMPLE_EXPECTED_RNOISE = {
    "data": np.array(
        [[0.52631587, 0.52631587], [0.23026317, 0.7236843]], dtype=np.float32
    ),
    "err": np.array([[14.712976, 14.712976], [14.711851, 14.713726]], dtype=np.float32),
    "var_poisson": np.array(
        [[0.05886428, 0.05886428], [0.02575312, 0.08093839]], dtype=np.float32
    ),
    "var_rnoise": np.array(
        [[216.4128, 216.4128], [216.4128, 216.4128]], dtype=np.float32
    ),
}


# #####
# Tests
# #####
def test_bad_readpattern():
    """Ensure error is raised on bad readpattern"""
    ramp_model, gain_model, readnoise_model, dark_model = make_data(
        SIMPLE_RESULTANTS, 1, 0.01, False
    )
    bad_pattern = ramp_model.meta.exposure.read_pattern.data[1:]
    ramp_model.meta.exposure.read_pattern = bad_pattern

    with pytest.raises(RuntimeError):
        RampFitStep.call(
            ramp_model,
            algorithm="ols_cas22",
            override_gain=gain_model,
            override_readnoise=readnoise_model,
        )


@pytest.mark.parametrize(
    "attribute",
    ["data", "err", "var_poisson", "var_rnoise"],
    ids=["data", "err", "var_poisson", "var_rnoise"],
)
def test_fits(fit_ramps, attribute):
    """Check slopes"""
    image_model, expected = fit_ramps

    value = getattr(image_model, attribute)
    np.testing.assert_allclose(value, expected[attribute], 1e-05)


# ########
# Fixtures
# ########
@pytest.fixture(
    scope="module",
    params=[
        pytest.param(
            (SIMPLE_RESULTANTS, 1, 0.01, False, SIMPLE_EXPECTED_DEFAULT, False),
            id="default",
        ),  # No gain or noise
        pytest.param(
            (SIMPLE_RESULTANTS, 5, 0.01, False, SIMPLE_EXPECTED_GAIN, False),
            id="extragain",
        ),  # Increase gain
        pytest.param(
            (SIMPLE_RESULTANTS, 1, 100.0, False, SIMPLE_EXPECTED_RNOISE, False),
            id="extranoise",
        ),  # Increase noise
        pytest.param(
            (SIMPLE_RESULTANTS, 1, 0.01, False, SIMPLE_EXPECTED_DEFAULT, True),
            id="default-jump",
        ),  # No gain or noise
        pytest.param(
            (SIMPLE_RESULTANTS, 5, 0.01, False, SIMPLE_EXPECTED_GAIN, True),
            id="extragain-jump",
        ),  # Increase gain
        pytest.param(
            (SIMPLE_RESULTANTS, 1, 100.0, False, SIMPLE_EXPECTED_RNOISE, True),
            id="extranoise-jump",
        ),  # Increase noise
    ],
)
def fit_ramps(request):
    """Test ramp fits"""
    resultants, ingain, rnoise, randomize, expected, use_jump = request.param
    ramp_model, gain_model, readnoise_model, dark_model = make_data(
        resultants, ingain, rnoise, randomize
    )

    out_model = RampFitStep.call(
        ramp_model,
        algorithm="ols_cas22",
        use_ramp_jump_detection=use_jump,
        override_gain=gain_model,
        override_readnoise=readnoise_model,
    )

    return out_model, expected
