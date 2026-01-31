"""Ramp Fitting tests involving MultiAccum Tables"""

import numpy as np
import pytest

from romancal.ramp_fitting import RampFitStep

from .common import SIMPLE_RESULTANTS, make_data

SIMPLE_EXPECTED_DEFAULT = {
    "data": np.array(
        [[0.5482434, 0.5482434], [0.21930099, 0.6579002]],
        dtype=np.float32,
    ),
    "err": np.array(
        [[0.24518493, 0.24518453], [0.1550718, 0.26858887]],
        dtype=np.float16,
    ),
    "var_poisson": np.array(
        [[0.06011445, 0.06011425], [0.02404606, 0.07213879]],
        dtype=np.float16,
    ),
    "var_rnoise": np.array(
        [[1.2022689e-06, 1.2022689e-06], [1.2022329e-06, 1.2022729e-06]],
        dtype=np.float16,
    ),
    "chisq": np.array([[0.4, 1.6], [1, 3]], dtype=np.float16),
}
SIMPLE_EXPECTED_GAIN = {
    "data": np.array(
        [[0.54823464, 0.54823464], [0.21931194, 0.32894737]], dtype=np.float32
    ),
    "err": np.array(
        [[0.10965369, 0.10965278], [0.0693583, 0.1040483]], dtype=np.float16
    ),
    "var_poisson": np.array(
        [[0.01202273, 0.01202253], [0.00480937, 0.01082064]], dtype=np.float16
    ),
    "var_rnoise": np.array(
        [[1.2021728e-06, 1.2021728e-06], [1.2019930e-06, 5.4103184e-06]],
        dtype=np.float16,
    ),
    "chisq": np.array([[2, 8], [5, 0]], dtype=np.float16),
}
SIMPLE_EXPECTED_RNOISE = {
    "data": np.array(
        [[0.5263179, 0.5263179], [0.2302627, 0.72367555]], dtype=np.float32
    ),
    "err": np.array([[10.405058, 10.405058], [10.403467, 10.406118]], dtype=np.float16),
    "var_poisson": np.array(
        [[0.05886434, 0.05886419], [0.025753, 0.0809375]], dtype=np.float16
    ),
    "var_rnoise": np.array(
        [[108.20637, 108.20637], [108.20637, 108.20637]], dtype=np.float16
    ),
    "chisq": np.array([[4e-5, 2.4e-4], [6e-5, 3.6e-4]], dtype=np.float16),
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

    with pytest.raises(ValueError):
        RampFitStep.call(
            ramp_model,
            algorithm="likely",
            override_gain=gain_model,
            override_readnoise=readnoise_model,
        )


@pytest.mark.parametrize(
    "attribute",
    ["data", "err", "var_poisson", "var_rnoise", "chisq"],
    ids=["data", "err", "var_poisson", "var_rnoise", "chisq"],
)
def test_fits(fit_ramps, attribute):
    """Check slopes"""
    image_model, expected = fit_ramps

    value = getattr(image_model, attribute)
    expected_value = expected[attribute]
    precision = 1e-5 if expected_value.dtype == np.float32 else 1e-3
    np.testing.assert_allclose(value, expected_value, precision)


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
        algorithm="likely",
        use_ramp_jump_detection=use_jump,
        override_gain=gain_model,
        override_readnoise=readnoise_model,
        include_var_rnoise=True,
    )

    return out_model, expected
