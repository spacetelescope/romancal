"""
Unit tests for the Roman source detection step code
"""

from copy import deepcopy

import crds
import numpy as np
import pytest
import roman_datamodels as rdm
from astropy import units as u
from astropy.modeling.models import Gaussian2D
from astropy.stats import mad_std
from astropy.table import QTable
from photutils.datasets import make_model_image
from photutils.psf import PSFPhotometry
from roman_datamodels.datamodels import ImageModel

from romancal.source_catalog.psf import (
    azimuthally_smooth,
    fit_psf_to_image_model,
    get_gridded_psf_model,
)

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
    mod = ImageModel.create_fake_data(shape=shape)
    mod.data = np.ones(shape, dtype=np.float32)
    mod.meta.filename = "filename"
    mod.meta.instrument["optical_element"] = "F087"

    # add noise to data
    if noise is not None:
        setup_rng = np.random.default_rng(seed or 19)
        mod.data = setup_rng.normal(scale=noise, size=shape).astype("float32")
        mod.err = noise * (np.ones(shape, dtype=np.float32) * u.DN / u.s).value

    # add dq array
    mod.dq = np.zeros(shape, dtype=np.uint32)

    crds_parameters = mod.get_crds_parameters()
    crds_ref_file = crds.getreferences(
        crds_parameters, reftypes=["epsf"], observatory="roman"
    )
    psf_ref_file = crds_ref_file["epsf"]
    mod["psf_ref_model"] = rdm.open(psf_ref_file)

    # compute gridded PSF model:
    psf_model = get_gridded_psf_model(mod["psf_ref_model"])

    return mod, psf_model


def add_sources(image_model, psf_model, x_true, y_true, flux_true, background=10):
    params_table = QTable()
    params_table["x_0"] = np.atleast_1d(x_true)
    params_table["y_0"] = np.atleast_1d(y_true)
    params_table["flux"] = np.atleast_1d(flux_true)

    shape = image_model.data.shape
    image = make_model_image(shape, psf_model, params_table, model_shape=(19, 19))
    image += rng.normal(background, image_model.err, size=shape)

    image_model.data = image * np.ones_like(image_model.data)


class TestPSFFitting:
    def setup_method(self):
        self.image_model, self.psf_model = setup_inputs(
            shape=image_model_shape,
        )

    @pytest.mark.parametrize(
        "dx, dy, true_flux",
        zip(
            rng.uniform(-1, 1, n_trials),
            rng.uniform(-1, 1, n_trials),
            np.geomspace(1_000, 100_000, n_trials),
            strict=False,
        ),
    )
    def test_psf_fit(self, dx, dy, true_flux):
        # generate an ImageModel
        image_model = deepcopy(self.image_model)

        # add synthetic sources to the ImageModel:
        true_x = image_model_shape[0] / 2 + dx
        true_y = image_model_shape[1] / 2 + dy
        add_sources(image_model, self.psf_model, true_x, true_y, true_flux)
        init_data_stddev = np.std(image_model.data)

        # fit the PSF to the ImageModel:
        results_table, photometry = fit_psf_to_image_model(
            image_model=image_model,
            photometry_cls=PSFPhotometry,
            psf_model=self.psf_model,
            x_init=true_x,
            y_init=true_y,
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
        approx_snr = true_flux / init_data_stddev
        approx_fwhm = 1
        approx_centroid_err = approx_fwhm / approx_snr

        # centroid err heuristic above is an underestimate, so we scale it up:
        scale_factor_approx = 2

        assert np.all(
            results_table["x_err"] < scale_factor_approx * approx_centroid_err
        )
        assert np.all(
            results_table["y_err"] < scale_factor_approx * approx_centroid_err
        )


def test_azimuthally_smooth():
    """Test azimuthally smoothing"""
    grid_x, grid_y = np.mgrid[0:201, 0:201]
    gauss_model = Gaussian2D(1.0, 100, 100, 20, 20)
    gauss = gauss_model(grid_x, grid_y)
    smoothed = azimuthally_smooth(gauss, oversample=1)
    delta = gauss - smoothed

    assert np.mean(delta) < 1.0e-6
    assert mad_std(delta) < 1.0e-8
