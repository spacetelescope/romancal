"""
 Unit tests for the Roman source detection step code
"""

from copy import deepcopy

import numpy as np
import pytest
from astropy import units as u

from romancal.lib.basic_utils import recarray_to_ndarray
from romancal.lib.tests.test_psf import add_sources, setup_inputs
from romancal.source_detection import SourceDetectionStep

n_sources = 10
image_model_shape = (100, 100)


class TestSourceDetection:
    def setup_method(self):
        self.image_model, self.webbpsf_config, self.psf_model = setup_inputs(
            shape=image_model_shape,
        )

    @pytest.mark.webbpsf
    def test_dao_vs_psf(self, seed=0):
        rng = np.random.default_rng(seed)
        image_model = deepcopy(self.image_model)

        n_models = 10
        amp_true = rng.normal(1e3, 100, n_models)
        along_diag = np.arange(100, 950, 90) / 1000 * image_model_shape[0]
        x_true = along_diag + rng.normal(scale=0.5, size=n_models)
        y_true = along_diag + rng.normal(scale=0.5, size=n_models)

        add_sources(image_model, self.psf_model, x_true, y_true, amp_true)

        source_detect = SourceDetectionStep()
        source_detect.scalar_threshold = 100
        source_detect.peakmax = None
        dao_result = source_detect.process(image_model)
        idx, x_dao, y_dao, amp_dao = recarray_to_ndarray(
            dao_result.meta.source_detection.tweakreg_catalog
        ).T

        # check that all injected targets are found by DAO:
        assert len(x_dao) == len(x_true)

        source_detect.fit_psf = True
        psf_result = source_detect.process(image_model)
        psf_catalog = psf_result.meta.source_detection.psf_catalog

        extract_columns = ["xcentroid", "x_err", "ycentroid", "y_err", "flux_fit"]
        x_psf, x_err, y_psf, y_err, amp_psf = psf_catalog[extract_columns].itercols()

        # check that typical PSF centroids are more accurate than DAO centroids:
        assert np.median(np.abs(x_dao - x_true)) > np.median(np.abs(x_psf - x_true))

        # check that the typical/worst PSF centroid is still within some tolerance:
        pixel_scale = 0.11 * u.arcsec / u.pix
        centroid_residuals = np.abs(x_psf - x_true) * u.pix * pixel_scale
        assert np.max(centroid_residuals) < 11 * u.mas
        assert np.median(centroid_residuals) < 3 * u.mas

        # check that typical residuals are consistent with their errors:
        assert np.median(np.abs(x_psf - x_true) / x_err) < 2
