"""
Unit tests for PSF matching functionality.
"""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.modeling.fitting import TRFLSQFitter
from astropy.modeling.models import Gaussian2D
from roman_datamodels.datamodels import ImageModel

from romancal.source_catalog.psf_matching import (
    create_psf_matched_image,
    get_filter_wavelength,
    get_reddest_filter,
)


@pytest.fixture
def mock_image_model():
    """
    Create a mock ImageModel for testing.
    """
    # Create a simple test image with a Gaussian-like source
    ny, nx = 100, 100
    y, x = np.mgrid[:ny, :nx]
    cy, cx = ny // 2, nx // 2

    # Create a simple Gaussian-like source
    sigma = 3.0
    gauss = Gaussian2D(
        amplitude=1.0, x_mean=cx, y_mean=cy, x_stddev=sigma, y_stddev=sigma
    )
    data = gauss(x, y)
    rng = np.random.default_rng(42)
    data += 0.001 * rng.standard_normal((ny, nx))  # Add small noise

    err = 0.01 * np.ones_like(data)

    # Use create_fake_data to get proper metadata structure
    model = ImageModel.create_fake_data(shape=(ny, nx))
    model.data[:] = data
    model.err[:] = err
    model.meta.instrument.optical_element = "F087"

    return model


@pytest.fixture
def mock_psf_ref_model():
    """
    Create a mock PSF reference model factory for testing.
    """

    def _create_mock_psf(filter_name, oversampling=4):
        """
        Create a mock PSF reference model.
        """

        class MockPSFRef:
            """
            Mock PSF reference model.
            """

            def __init__(self, filter_name, oversampling):
                # Set metadata
                self.meta = type("Meta", (), {})()
                self.meta.instrument = type("Instrument", (), {})()
                self.meta.instrument.optical_element = filter_name
                self.oversample = oversampling
                self.oversampling = oversampling

                # Create a simple PSF stamp for create_l3_psf_model
                stamp_size = 41
                y, x = np.mgrid[:stamp_size, :stamp_size]
                cy, cx = stamp_size // 2, stamp_size // 2

                # Create Gaussian PSF - broader for longer wavelengths
                wavelength = get_filter_wavelength(filter_name)
                # Scale FWHM with wavelength
                fwhm = 2.0 + wavelength  # Simple scaling
                sigma = fwhm / 2.355
                psf = np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma**2))
                psf = psf / psf.sum()

                self.psf_data = psf

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

        return MockPSFRef(filter_name, oversampling)

    return _create_mock_psf


@pytest.fixture
def mock_create_l3_psf_patch():
    """
    Create a mock patch for create_l3_psf_model.

    Returns a context manager that can be used to patch create_l3_psf_model
    in tests that need to perform PSF matching.
    """

    def mock_create_l3_psf(psf_ref_model):
        """Mock create_l3_psf_model."""

        class MockL3PSF:
            def __init__(self):
                self.data = psf_ref_model.psf_data
                self.oversampling = np.array([4, 4])

        return MockL3PSF()

    return patch(
        "romancal.source_catalog.psf_matching.create_l3_psf_model",
        side_effect=mock_create_l3_psf,
    )


def test_get_filter_wavelength():
    """
    Test filter wavelength extraction.
    """
    assert get_filter_wavelength("F062") == 0.62
    assert get_filter_wavelength("F087") == 0.87
    assert get_filter_wavelength("F106") == 1.06
    assert get_filter_wavelength("F129") == 1.29
    assert get_filter_wavelength("F158") == 1.58
    assert get_filter_wavelength("F184") == 1.84
    assert get_filter_wavelength("F213") == 2.13

    # Test with 'm' suffix (PSF-matched column names)
    assert get_filter_wavelength("F158m") == 1.58
    assert get_filter_wavelength("f184m") == 1.84

    # Test invalid filter names
    assert get_filter_wavelength("invalid") == 0
    assert get_filter_wavelength("") == 0


def test_create_psf_matched_image_basic(
    mock_image_model, mock_psf_ref_model, mock_create_l3_psf_patch
):
    """
    Test basic PSF matching functionality.
    """
    # Create input and target PSF models
    input_psf_ref = mock_psf_ref_model("F087")  # Narrower PSF
    target_psf_ref = mock_psf_ref_model("F184")  # Broader PSF

    with mock_create_l3_psf_patch:
        # Create matched image
        matched_model = create_psf_matched_image(
            mock_image_model, input_psf_ref, target_psf_ref
        )

    # Check that output is an ImageModel
    assert isinstance(matched_model, ImageModel)

    # Check that data shape is preserved
    assert matched_model.data.shape == mock_image_model.data.shape

    # Check that errors are propagated
    assert matched_model.err is not None
    assert matched_model.err.shape == mock_image_model.err.shape
    assert np.all(matched_model.err > 0)
    assert np.all(np.isfinite(matched_model.err))


def test_create_psf_matched_image_same_filter(mock_image_model, mock_psf_ref_model):
    """
    Test that PSF matching is skipped when filters are the same.
    """
    # Create PSF models with same filter
    input_psf_ref = mock_psf_ref_model("F087")
    target_psf_ref = mock_psf_ref_model("F087")

    matched_model = create_psf_matched_image(
        mock_image_model, input_psf_ref, target_psf_ref
    )

    # Check that matching was skipped - returns the input model
    assert matched_model is mock_image_model


def test_create_psf_matched_image_broader_input(mock_image_model, mock_psf_ref_model):
    """
    Test that PSF matching is skipped when input PSF is already broader.
    """
    # Create PSF models where input is broader than target
    input_psf_ref = mock_psf_ref_model("F184")  # Broader
    target_psf_ref = mock_psf_ref_model("F087")  # Narrower

    mock_image_model.meta.instrument.optical_element = "F184"

    matched_model = create_psf_matched_image(
        mock_image_model, input_psf_ref, target_psf_ref
    )

    # Check that matching was skipped - returns the input model
    assert matched_model is mock_image_model


def test_create_psf_matched_image_params(mock_image_model, mock_psf_ref_model):
    """
    Test that all required parameters must be provided.
    """
    input_psf_ref = mock_psf_ref_model("F087")

    match = "missing 1 required positional argument"
    with pytest.raises(TypeError, match=match):
        create_psf_matched_image(mock_image_model, input_psf_ref)


def test_create_psf_matched_image_flux_conservation(
    mock_image_model, mock_psf_ref_model, mock_create_l3_psf_patch
):
    """
    Test that total flux is roughly preserved after PSF matching.
    """
    input_psf_ref = mock_psf_ref_model("F087")
    target_psf_ref = mock_psf_ref_model("F184")

    input_flux = np.sum(mock_image_model.data)

    with mock_create_l3_psf_patch:
        matched_model = create_psf_matched_image(
            mock_image_model, input_psf_ref, target_psf_ref
        )

    matched_flux = np.sum(matched_model.data)

    # Flux should be approximately preserved
    assert np.isclose(matched_flux, input_flux, rtol=0.05)


def test_get_reddest_filter():
    """
    Test automatic reddest filter selection.
    """

    # Create a mock library with multiple filters
    class MockLibrary:
        def __init__(self, models):
            self.models = models
            self._index = 0

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def __iter__(self):
            return iter(self.models)

        def shelve(self, model, modify=False):
            pass

    # Create models with different filters
    model_f087 = ImageModel.create_fake_data(shape=(10, 10))
    model_f087.meta.instrument.optical_element = "F087"

    model_f129 = ImageModel.create_fake_data(shape=(10, 10))
    model_f129.meta.instrument.optical_element = "F129"

    model_f184 = ImageModel.create_fake_data(shape=(10, 10))
    model_f184.meta.instrument.optical_element = "F184"

    library = MockLibrary([model_f087, model_f129, model_f184])

    reddest_filter = get_reddest_filter(library)

    # Should select F184 (longest wavelength)
    assert reddest_filter == "F184"


def test_create_psf_matched_image_invalid_input(mock_psf_ref_model):
    """
    Test that invalid model type raises ValueError.
    """
    invalid_model = "not a model"
    input_psf_ref = mock_psf_ref_model("F087")
    target_psf_ref = mock_psf_ref_model("F184")

    match = "model must be an ImageModel or MosaicModel"
    with pytest.raises(ValueError, match=match):
        create_psf_matched_image(invalid_model, input_psf_ref, target_psf_ref)


def test_psf_matching_kernel_validation():
    """
    Validation test for PSF matching with known Gaussian PSFs.

    This test creates a noiseless image with a single centered Gaussian
    source (sigma=3) and PSF (sigma=3). The target PSF is a Gaussian
    with sigma=5. The matching kernel should be a Gaussian with sigma=4,
    and the final PSF-matched source should have sigma=5.

    This validates that:
    - The matching kernel has the expected sigma=4
    - The output source has the expected sigma=5
    """
    # Create noiseless image with centered Gaussian source
    ny, nx = 201, 201  # Larger to avoid edge effects
    y, x = np.mgrid[:ny, :nx]
    cy, cx = ny // 2, nx // 2

    # Source with sigma=3
    source_sigma = 3.0
    gauss_source = Gaussian2D(
        amplitude=1000.0,
        x_mean=cx,
        y_mean=cy,
        x_stddev=source_sigma,
        y_stddev=source_sigma,
    )
    data = gauss_source(x, y)

    # Create model
    model = ImageModel.create_fake_data(shape=(ny, nx))
    model.data[:] = data
    model.err = 0.01 * np.ones_like(data)
    model.meta.instrument.optical_element = "F087"

    # Create PSF models with known Gaussian PSFs
    # Note: PSFs will be oversampled by a factor of 4
    oversample = 4

    class MockPSFRef:
        """Mock PSF reference model with Gaussian PSF."""

        def __init__(self, filter_name, sigma):
            self.meta = type("Meta", (), {})()
            self.meta.instrument = type("Instrument", (), {})()
            self.meta.instrument.optical_element = filter_name
            self.oversample = oversample
            self.oversampling = oversample

            # Create Gaussian PSF stamp (oversampled)
            # sigma in the oversampled space
            sigma_oversampled = sigma * oversample
            stamp_size = 41 * oversample  # Larger stamp for oversampled PSF
            y_psf, x_psf = np.mgrid[:stamp_size, :stamp_size]
            cy_psf, cx_psf = stamp_size // 2, stamp_size // 2

            gauss_psf = Gaussian2D(
                amplitude=1.0,
                x_mean=cx_psf,
                y_mean=cy_psf,
                x_stddev=sigma_oversampled,
                y_stddev=sigma_oversampled,
            )
            psf = gauss_psf(x_psf, y_psf)
            psf = psf / psf.sum()  # Normalize

            self.psf_data = psf

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    # Input PSF: sigma=3 (same as source)
    input_psf_sigma = 3.0
    input_psf_ref = MockPSFRef("F087", input_psf_sigma)

    # Target PSF: sigma=5
    target_psf_sigma = 5.0
    target_psf_ref = MockPSFRef("F184", target_psf_sigma)

    # Expected matching kernel sigma: sqrt(5^2 - 3^2) = sqrt(25 - 9) = 4
    expected_kernel_sigma = np.sqrt(target_psf_sigma**2 - input_psf_sigma**2)

    # Mock create_l3_psf_model to capture the kernel
    captured_kernel = {}

    def mock_create_l3_psf(psf_ref_model):
        """Mock create_l3_psf_model."""

        class MockL3PSF:
            def __init__(self):
                self.data = psf_ref_model.psf_data
                self.oversampling = np.array([4, 4])

        return MockL3PSF()

    # Also patch create_convolution_kernel to capture it
    def mock_create_kernel(input_psf, target_psf, downsample=1):
        """Mock that captures the kernel."""
        # Import the real function
        from romancal.source_catalog.psf import create_convolution_kernel

        kernel = create_convolution_kernel(input_psf, target_psf, downsample=downsample)
        captured_kernel["kernel"] = kernel
        return kernel

    with (
        patch(
            "romancal.source_catalog.psf_matching.create_l3_psf_model",
            side_effect=mock_create_l3_psf,
        ),
        patch(
            "romancal.source_catalog.psf_matching.create_convolution_kernel",
            side_effect=mock_create_kernel,
        ),
    ):
        matched_model = create_psf_matched_image(model, input_psf_ref, target_psf_ref)

    # Check that we captured the kernel
    assert "kernel" in captured_kernel
    kernel = captured_kernel["kernel"]

    # Fit a 2D Gaussian to the matching kernel to extract sigma
    y_k, x_k = np.mgrid[: kernel.shape[0], : kernel.shape[1]]
    cy_k, cx_k = kernel.shape[0] // 2, kernel.shape[1] // 2

    # Initial guess for Gaussian fit
    gauss_init = Gaussian2D(
        amplitude=kernel.max(),
        x_mean=cx_k,
        y_mean=cy_k,
        x_stddev=expected_kernel_sigma,
        y_stddev=expected_kernel_sigma,
    )

    fitter = TRFLSQFitter()
    gauss_fit = fitter(gauss_init, x_k, y_k, kernel)

    # Check kernel sigma (should be ~4)
    fitted_kernel_sigma = (gauss_fit.x_stddev.value + gauss_fit.y_stddev.value) / 2
    assert np.isclose(fitted_kernel_sigma, expected_kernel_sigma, rtol=0.01), (
        f"Kernel sigma {fitted_kernel_sigma:.2f} != "
        f"expected {expected_kernel_sigma:.2f}"
    )

    # Fit a 2D Gaussian to the output source to extract sigma
    gauss_output_init = Gaussian2D(
        amplitude=matched_model.data.max(),
        x_mean=cx,
        y_mean=cy,
        x_stddev=target_psf_sigma,
        y_stddev=target_psf_sigma,
    )

    gauss_output_fit = fitter(gauss_output_init, x, y, matched_model.data)

    # Check output source sigma (should be ~5)
    fitted_output_sigma = (
        gauss_output_fit.x_stddev.value + gauss_output_fit.y_stddev.value
    ) / 2
    assert np.isclose(fitted_output_sigma, target_psf_sigma, rtol=0.01), (
        f"Output source sigma {fitted_output_sigma:.2f} != "
        f"expected {target_psf_sigma:.2f}"
    )
