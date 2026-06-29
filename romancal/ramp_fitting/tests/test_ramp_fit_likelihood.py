"""Ramp Fitting tests involving MultiAccum Tables"""

import numpy as np
import pytest
from roman_datamodels.dqflags import pixel

from romancal.ramp_fitting import RampFitStep

from .common import (
    ROMAN_READ_TIME,
    SIMPLE_RESULTANTS,
    create_linear_ramp,
    generate_wfi_reffiles,
    make_data,
    model_from_resultants,
)

SIMPLE_EXPECTED_DEFAULT = {
    "data": np.array(
        [[0.54824126, 0.54824126], [0.21930373, 0.6579057]],
        dtype=np.float32,
    ),
    "err": np.array(
        [[0.2452, 0.2452], [0.155, 0.2686]],
        dtype=np.float16,
    ),
    "var_poisson": np.array(
        [[0.06012, 0.06012], [0.02405, 0.07214]],
        dtype=np.float16,
    ),
    "var_rnoise": np.array(
        [[2.4e-06, 2.4e-06], [2.4e-06, 2.4e-06]],
        dtype=np.float16,
    ),
    "chisq": np.array([[0.4, 1.6], [0.9995, 3.0]], dtype=np.float16),
}
SIMPLE_EXPECTED_GAIN = {
    "data": np.array(
        [[0.54822373, 0.54822373], [0.21932559, 0.32894737]], dtype=np.float32
    ),
    "err": np.array([[0.1097, 0.1097], [0.0694, 0.10406]], dtype=np.float16),
    "var_poisson": np.array(
        [[0.012024, 0.012024], [0.00481, 0.01082]], dtype=np.float16
    ),
    "var_rnoise": np.array(
        [[2.384e-06, 2.384e-06], [2.384e-06, 1.085e-05]],
        dtype=np.float16,
    ),
    "chisq": np.array([[1.998, 7.996], [4.992, 0.0]], dtype=np.float16),
}
SIMPLE_EXPECTED_RNOISE = {
    "data": np.array(
        [[0.5263168, 0.5263168], [0.23026292, 0.72367984]], dtype=np.float32
    ),
    "err": np.array([[14.71, 14.71], [14.71, 14.71]], dtype=np.float16),
    "var_poisson": np.array([[0.05887, 0.05887], [0.02576, 0.08093]], dtype=np.float16),
    "var_rnoise": np.array([[216.4, 216.4], [216.4, 216.4]], dtype=np.float16),
    "chisq": np.array([[2.0e-05, 1.2e-04], [3.0e-05, 1.8e-04]], dtype=np.float16),
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


def test_flag_large_events_withsnowball():
    """Test that large events are flagged"""
    resultants = create_linear_ramp(n_resultants=20, nrows=100, ncols=100)
    model, m_gain, m_rnoise, m_dark = make_data(resultants, 6, 0.01, False)

    # square of saturation surrounded by jump -> snowball
    # 112 pixels (121 minus 9) initially have a jump.
    model.data[10:, 46:57, 46:57] += 300
    model.data[10:, 50:53, 50:53] = 1e5
    model.groupdq[10:, 50:53, 50:53] = pixel.SATURATED

    m_image = RampFitStep.call(
        model,
        algorithm="likely",
        override_gain=m_gain,
        override_readnoise=m_rnoise,
        expand_large_events=True,
        min_sat_area=1,
        min_jump_area=6,
        expand_factor=1.9,
        edge_size=0,
        sat_required_snowball=True,
        min_sat_radius_extend=0.5,
        sat_expand=2,
    )
    assert np.std(m_image.data) < 1e-5
    n_jump_expanded = np.sum(m_image.dq == pixel.JUMP_DET)

    # Check that the uncertainties are the same for all pixels with jumps
    # and save this average uncertainty
    meanerr_jumppixels_new = np.mean(m_image.err[m_image.dq == pixel.JUMP_DET])
    assert np.std(m_image.err[m_image.dq == pixel.JUMP_DET]) < 1e-6

    # Now without snowball flagging
    model, m_gain, m_rnoise, m_dark = make_data(resultants, 6, 0.01, False)
    model.data[10:, 46:57, 46:57] += 300
    model.data[10:, 50:53, 50:53] = 1e5
    model.groupdq[10:, 50:53, 50:53] = pixel.SATURATED

    m_image = RampFitStep.call(
        model,
        algorithm="likely",
        override_gain=m_gain,
        override_readnoise=m_rnoise,
        expand_large_events=False,
    )
    assert np.std(m_image.data) < 1e-5
    n_jump_original = np.sum(m_image.dq == pixel.JUMP_DET)

    # Check that the uncertainties are the same for all pixels with jumps
    # originally flagged.  Check that this uncertainty matches the new value.

    meanerr_jumppixels_orig = np.mean(m_image.err[m_image.dq == pixel.JUMP_DET])
    assert np.std(m_image.err[m_image.dq == pixel.JUMP_DET]) < 1e-6
    assert np.abs(meanerr_jumppixels_orig - meanerr_jumppixels_new) < 1e-5

    assert n_jump_original == 112 and n_jump_expanded > 300 and n_jump_expanded < 600


def test_record_jumps():
    """jump_indices records which resultant each jump was flagged in."""
    resultants = create_linear_ramp(n_resultants=20, nrows=20, ncols=20)
    model, m_gain, m_rnoise, m_dark = make_data(resultants, 6, 0.01, False)

    # Inject a clear jump at resultant 10 in a single science pixel.
    jump_group, row, col = 10, 8, 12
    model.data[jump_group:, row, col] += 300

    m_image = RampFitStep.call(
        model,
        algorithm="likely",
        override_gain=m_gain,
        override_readnoise=m_rnoise,
        record_jumps=True,
    )

    # The pixel is flagged JUMP_DET in the trimmed 2D dq.
    assert m_image.dq[row - 4, col - 4] & pixel.JUMP_DET

    jidx = m_image.jump_indices
    assert jidx.ndim == 2 and jidx.shape[1] == 3
    # The injected jump is recorded at the trimmed-frame pixel and resultant 10.
    match = (jidx[:, 1] == row - 4) & (jidx[:, 2] == col - 4)
    assert match.any()
    assert jump_group in jidx[match, 0]


def test_record_jumps_disabled():
    """record_jumps=False omits the jump_indices extension."""
    resultants = create_linear_ramp(n_resultants=20, nrows=20, ncols=20)
    model, m_gain, m_rnoise, m_dark = make_data(resultants, 6, 0.01, False)

    m_image = RampFitStep.call(
        model,
        algorithm="likely",
        override_gain=m_gain,
        override_readnoise=m_rnoise,
        record_jumps=False,
    )

    assert "jump_indices" not in m_image


def test_uniformweighting():
    """Ensure uniform weighting slope only matches in the read noise limit"""
    ramp_model, gain_model, readnoise_model, dark_model = make_data(
        SIMPLE_RESULTANTS, 0.01, 1000, False
    )

    out_model = RampFitStep.call(
        ramp_model,
        algorithm="likely",
        override_gain=gain_model,
        override_readnoise=readnoise_model,
        include_var_rnoise=True,
    )

    # Fully in the read noise limit: the uniform-weighting and
    # optimal-weighting sloeps should be the same.

    np.testing.assert_allclose(out_model.dumo, 0, atol=1e-5)

    ramp_model, gain_model, readnoise_model, dark_model = make_data(
        SIMPLE_RESULTANTS, 10, 0.01, False
    )

    out_model = RampFitStep.call(
        ramp_model,
        algorithm="likely",
        override_gain=gain_model,
        override_readnoise=readnoise_model,
        include_var_rnoise=True,
    )

    # Now we are in the photon noise limit.  The uniform-weighting and
    # optimal-weighting slopes should be different.

    with pytest.raises(AssertionError):
        np.testing.assert_allclose(out_model.dumo, 0, atol=1e-2)


@pytest.mark.parametrize("algorithm", ["ols_cas22", "likely"])
def test_uncertainty_matches_scatter(algorithm):
    """Verify that uncertainties match the observed scatter in slopes."""
    rng = np.random.default_rng(42)
    ntrial, readnoise, flux = 2000, 20.0, 10.0
    # Non-trivial read pattern
    read_pattern = [[1, 2, 3, 4], [5], [6, 7, 8], [9, 10, 11, 12]]
    nreads = np.array([len(g) for g in read_pattern])
    n_resultants = len(read_pattern)

    # Generate counts in each read
    total_reads = sum(nreads)
    per_read = rng.poisson(flux, size=(total_reads, ntrial)).astype(np.float32)
    cumulative = np.cumsum(per_read, axis=0)
    # Average all read values within each resultant group
    boundaries = np.concatenate([[0], np.cumsum(nreads)])
    resultants = np.array(
        [
            np.mean(cumulative[boundaries[i] : boundaries[i + 1]], axis=0)
            for i in range(n_resultants)
        ],
        dtype=np.float32,
    )
    # Read noise averages down as sqrt(n) within each resultant
    resultants += (
        rng.normal(0, 1, size=(n_resultants, ntrial)).astype(np.float32)
        * (readnoise / np.sqrt(nreads))[:, np.newaxis]
    )

    resultants_3d = resultants[:, np.newaxis, :]
    ramp_model = model_from_resultants(resultants_3d, read_pattern=read_pattern)
    gain_model, readnoise_model, _ = generate_wfi_reffiles(
        ramp_model.shape[1:], ingain=1, rnoise=readnoise, randomize=False
    )

    out = RampFitStep.call(
        ramp_model,
        algorithm=algorithm,
        override_gain=gain_model,
        override_readnoise=readnoise_model,
    )

    # out.data shape is [1, ntrial]: one science row, ntrial columns
    slopes = out.data[0]
    true_rate = flux / ROMAN_READ_TIME  # simulation flux is electrons/frame
    residuals = (slopes - true_rate) / out.err[0]
    np.testing.assert_allclose(np.std(slopes), np.mean(out.err[0]), rtol=0.05)
    np.testing.assert_allclose(np.mean(residuals), 0, atol=0.1)
    np.testing.assert_allclose(np.std(residuals), 1.0, rtol=0.05)


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
