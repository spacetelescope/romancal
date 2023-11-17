"""
 Unit tests for the Roman source detection step code
"""

import os
import tempfile
from copy import deepcopy

import numpy as np
import pytest
from astropy import units as u
from photutils.psf import PSFPhotometry
from roman_datamodels import maker_utils as testutil
from roman_datamodels.datamodels import ImageModel

from romancal.lib.psf import create_gridded_psf_model, fit_psf_to_image_model

n_trials = 15
image_model_shape = (50, 50)
rng = np.random.default_rng(0)


def setup_inputs(
    shape=image_model_shape,
    noise=1.0,
    seed=None,
):
    """
    Return ImageModel of level 2 image.
    """
    wfi_image = testutil.mk_level2_image(shape=shape)
    wfi_image.data = u.Quantity(
        np.ones(shape, dtype=np.float32), u.electron / u.s, dtype=np.float32
    )
    wfi_image.meta.filename = "filename"
    wfi_image.meta.instrument["optical_element"] = "F087"

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

    filt = mod.meta.instrument["optical_element"]
    detector = mod.meta["instrument"]["detector"].replace("WFI", "SCA")

    # input parameters for PSF model:
    webbpsf_config = dict(
        filt=filt,
        detector=detector,
        oversample=12,
        fov_pixels=15,
    )

    dir_path = tempfile.gettempdir()
    filename_prefix = f"psf_model_{webbpsf_config['filt']}"
    file_path = os.path.join(dir_path, filename_prefix)

    # compute gridded PSF model:
    psf_model, centroids = create_gridded_psf_model(
        file_path,
        webbpsf_config["filt"],
        webbpsf_config["detector"],
        oversample=webbpsf_config["oversample"],
        fov_pixels=webbpsf_config["fov_pixels"],
        overwrite=True,
        logging_level="ERROR",
    )

    return mod, webbpsf_config, psf_model


def add_sources(image_model, psf_model, x_true, y_true, amp_true, background=10):
    shape = image_model.data.shape

    x_true, y_true, amp_true = (
        np.atleast_1d(val) for val in [x_true, y_true, amp_true]
    )

    fit_models = []
    for amp, x, y in zip(amp_true, x_true, y_true):
        model = deepcopy(psf_model)
        model.parameters = [amp, x, y]
        fit_models.append(model)

    photometry = PSFPhotometry(psf_model, (15, 15))
    photometry.fit_results = dict(local_bkg=np.ones(len(x_true)) * 0)
    photometry._fit_models = fit_models

    image = photometry.make_model_image(shape, (19, 19))

    image += rng.normal(background, 1, size=shape)
    image_model.data = image * np.ones_like(image_model.data)
    image_model.err = background * np.ones_like(image_model.err)


class TestPSFFitting:
    def setup_method(self):
        self.image_model, self.webbpsf_config, self.psf_model = setup_inputs(
            shape=image_model_shape,
        )

    @pytest.mark.webbpsf
    @pytest.mark.parametrize(
        "dx, dy, true_amp",
        zip(
            rng.uniform(-1, 1, n_trials),
            rng.uniform(-1, 1, n_trials),
            np.geomspace(1_000, 100_000, n_trials),
        ),
    )
    def test_psf_fit(self, dx, dy, true_amp):
        # generate an ImageModel
        image_model = deepcopy(self.image_model)
        init_data_stddev = np.std(image_model.data.value)

        # add synthetic sources to the ImageModel:
        true_x = image_model_shape[0] / 2 + dx
        true_y = image_model_shape[1] / 2 + dy
        add_sources(image_model, self.psf_model, true_x, true_y, true_amp)

        if self.webbpsf_config["fov_pixels"] % 2 == 0:
            fit_shape = (
                self.webbpsf_config["fov_pixels"] + 1,
                self.webbpsf_config["fov_pixels"] + 1,
            )
        else:
            fit_shape = (
                self.webbpsf_config["fov_pixels"],
                self.webbpsf_config["fov_pixels"],
            )

        # fit the PSF to the ImageModel:
        results_table, photometry = fit_psf_to_image_model(
            image_model=image_model,
            photometry_cls=PSFPhotometry,
            psf_model=self.psf_model,
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
        # the above test to pass even when the fits are bad. Use overly-simple
        # approximation that astrometric uncertainty be proportional to the
        # PSF's FWHM / SNR:
        approx_snr = true_amp / init_data_stddev
        approx_fwhm = 1
        approx_centroid_err = approx_fwhm / approx_snr

        # centroid err heuristic above is an underestimate, so we scale it up:
        scale_factor_approx = 100

        assert np.all(
            results_table["x_err"] < scale_factor_approx * approx_centroid_err
        )
        assert np.all(
            results_table["y_err"] < scale_factor_approx * approx_centroid_err
        )
