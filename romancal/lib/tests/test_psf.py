"""
 Unit tests for the Roman source detection step code
"""

import os
import tempfile

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import overlap_slices
from photutils.psf import PSFPhotometry
from roman_datamodels import maker_utils as testutil
from roman_datamodels.datamodels import ImageModel

from romancal.lib.psf import create_gridded_psf_model, fit_psf_to_image_model

n_sources = 10
image_model_shape = (100, 100)
rng = np.random.default_rng(0)


@pytest.fixture
def setup_inputs():
    def _setup(
        nrows=image_model_shape[0], ncols=image_model_shape[1], noise=1.0, seed=None
    ):
        """
        Return ImageModel of level 2 image.
        """
        shape = (nrows, ncols)
        wfi_image = testutil.mk_level2_image(shape=shape)
        wfi_image.data = u.Quantity(
            np.ones(shape, dtype=np.float32), u.electron / u.s, dtype=np.float32
        )
        wfi_image.meta.filename = "filename"

        # add noise to data
        if noise is not None:
            setup_rng = np.random.default_rng(seed or 19)
            wfi_image.data = u.Quantity(
                setup_rng.normal(scale=noise, size=shape),
                u.electron / u.s,
                dtype=np.float32,
            )
            wfi_image.err = noise * np.ones(shape, dtype=np.float32) * u.electron / u.s

        # add dq array
        wfi_image.dq = np.zeros(shape, dtype=np.uint32)

        # construct ImageModel
        mod = ImageModel(wfi_image)

        return mod

    return _setup


def add_synthetic_sources(
    image_model,
    psf_model,
    true_x,
    true_y,
    true_amp,
    oversample,
    xname="x_0",
    yname="y_0",
):
    fit_models = []

    # ensure truths are arrays:
    true_x, true_y, true_amp = (
        np.atleast_1d(truth) for truth in [true_x, true_y, true_amp]
    )

    for x, y, amp in zip(true_x, true_y, true_amp):
        psf = psf_model.copy()
        psf.parameters = [amp, x, y]
        fit_models.append(psf)

    synth_image = image_model.data
    synth_err = image_model.err
    psf_shape = np.array(psf_model.data.shape[1:]) // oversample

    for fit_model in fit_models:
        x0 = getattr(fit_model, xname).value
        y0 = getattr(fit_model, yname).value
        slc_lg, _ = overlap_slices(synth_image.shape, psf_shape, (y0, x0), mode="trim")
        yy, xx = np.mgrid[slc_lg]
        model_data = fit_model(xx, yy) * image_model.data.unit
        model_err = np.sqrt(model_data.value) * model_data.unit
        synth_image[slc_lg] += (
            np.random.normal(
                model_data.to_value(image_model.data.unit),
                model_err.to_value(image_model.data.unit),
                size=model_data.shape,
            )
            * image_model.data.unit
        )
        synth_err[slc_lg] = np.sqrt(synth_err[slc_lg] ** 2 + model_err**2)


@pytest.mark.parametrize(
    "dx, dy, true_amp",
    zip(
        rng.uniform(-1, 1, n_sources),
        rng.uniform(-1, 1, n_sources),
        np.geomspace(10, 10_000, n_sources),
    ),
)
def test_psf_fit(setup_inputs, dx, dy, true_amp, seed=42):
    # input parameters for PSF model:
    filt = "F087"
    detector = "SCA01"
    oversample = 12
    fov_pixels = 15

    dir_path = tempfile.gettempdir()
    filename_prefix = f"psf_model_{filt}"
    file_path = os.path.join(dir_path, filename_prefix)

    # compute gridded PSF model:
    psf_model, centroids = create_gridded_psf_model(
        file_path,
        filt,
        detector,
        oversample=oversample,
        fov_pixels=fov_pixels,
        overwrite=False,
        logging_level="ERROR",
    )

    # generate an ImageModel
    image_model = setup_inputs(seed=seed)
    init_data_stddev = np.std(image_model.data.value)

    # add synthetic sources to the ImageModel:
    true_x = image_model_shape[0] / 2 + dx
    true_y = image_model_shape[1] / 2 + dy
    add_synthetic_sources(
        image_model, psf_model, true_x, true_y, true_amp, oversample=oversample
    )

    if fov_pixels % 2 == 0:
        fit_shape = (fov_pixels + 1, fov_pixels + 1)
    else:
        fit_shape = (fov_pixels, fov_pixels)

    # fit the PSF to the ImageModel:
    results_table, photometry = fit_psf_to_image_model(
        image_model=image_model,
        photometry_cls=PSFPhotometry,
        psf_model=psf_model,
        x_init=true_x,
        y_init=true_y,
        fit_shape=fit_shape,
    )

    # difference between input and output, normalized by the
    # uncertainty. Has units of sigma:
    delta_x = np.abs(true_x - results_table["x_fit"]) / results_table["x_err"]
    delta_y = np.abs(true_y - results_table["y_fit"]) / results_table["y_err"]

    sigma_threshold = 3.5
    assert np.all(delta_x < sigma_threshold)
    assert np.all(delta_y < sigma_threshold)

    # now check that the uncertainties aren't way too large, which could cause
    # the above test to pass even when the fits are bad. Use overly-simple approximation
    # that astrometric uncertainty be proportional to the PSF's FWHM / SNR:
    approx_snr = true_amp / init_data_stddev
    approx_fwhm = 1
    approx_centroid_err = approx_fwhm / approx_snr

    # centroid err heuristic above is an underestimate, so we scale it up:
    scale_factor_approx = 100

    assert np.all(results_table["x_err"] < scale_factor_approx * approx_centroid_err)
    assert np.all(results_table["y_err"] < scale_factor_approx * approx_centroid_err)
