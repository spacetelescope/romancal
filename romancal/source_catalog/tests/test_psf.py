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
from astropy.convolution import convolve

from romancal.source_catalog import psf
from romancal.source_catalog.psf import (
    fit_psf_to_image_model,
    get_gridded_psf_model,
)

n_trials = 15
image_model_shape = (50, 50)
rng = np.random.default_rng(0)


@pytest.fixture(scope="module")
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
    crds_parameters_2 = deepcopy(crds_parameters)
    crds_ref_file = crds.getreferences(
        crds_parameters, reftypes=["epsf"], observatory="roman"
    )
    psf_ref_file = crds_ref_file["epsf"]
    psf_ref_model_f087 = rdm.open(psf_ref_file)

    # compute gridded PSF model:
    psf_model = get_gridded_psf_model(psf_ref_model_f087)

    crds_parameters["meta.instrument.optical_element"] = "F184"
    crds_ref_file = crds.getreferences(
        crds_parameters, reftypes=["epsf"], observatory="roman"
    )
    psf_ref_file = crds_ref_file["epsf"]
    psf_ref_model_f184 = rdm.open(psf_ref_file)

    return dict(image=mod,
                psf_model=psf_model,
                psf_ref_model_f087=psf_ref_model_f087,
                psf_ref_model_f184=psf_ref_model_f184)


def add_sources(image_model, psf_model, x_true, y_true, flux_true, background=10):
    params_table = QTable()
    params_table["x_0"] = np.atleast_1d(x_true)
    params_table["y_0"] = np.atleast_1d(y_true)
    params_table["flux"] = np.atleast_1d(flux_true)

    shape = image_model.data.shape
    image = make_model_image(shape, psf_model, params_table, model_shape=(19, 19))
    image += rng.normal(background, image_model.err, size=shape)

    image_model.data = image * np.ones_like(image_model.data)


@pytest.mark.parametrize(
    "dx, dy, true_flux",
    zip(
        rng.uniform(-1, 1, n_trials),
        rng.uniform(-1, 1, n_trials),
        np.geomspace(1_000, 100_000, n_trials),
        strict=False,
    ),
)
def test_psf_fit(setup_inputs, dx, dy, true_flux):
    image_model = setup_inputs["image"]
    psf_model = setup_inputs["psf_model"]
    image_model = deepcopy(image_model)

    # add synthetic sources to the ImageModel:
    true_x = image_model_shape[0] / 2 + dx
    true_y = image_model_shape[1] / 2 + dy
    add_sources(image_model, psf_model, true_x, true_y, true_flux)
    init_data_stddev = np.std(image_model.data)

    # fit the PSF to the ImageModel:
    results_table, photometry = fit_psf_to_image_model(
        image_model=image_model,
        photometry_cls=PSFPhotometry,
        psf_model=psf_model,
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

    assert np.all(results_table["x_err"] < scale_factor_approx * approx_centroid_err)
    assert np.all(results_table["y_err"] < scale_factor_approx * approx_centroid_err)


# new routines: render_stamp, _get_jitter_params, _evaluate_gaussian_fft, add_jitter
def test_render_stamp(setup_inputs):
    # some basic tests that we can render a psf
    psf_model = setup_inputs["psf_model"]
    stamp = psf.render_stamp(2000, 2000, psf_model, 19)
    assert 0.9 < np.sum(stamp) < 1.1
    assert stamp.shape[0] == stamp.shape[1] == 19
    stamp = psf.render_stamp(0, 0, psf_model, 19)
    assert 0.9 < np.sum(stamp) < 1.1
    stamp = psf.render_stamp(0, 0, psf_model, 1)
    assert np.sum(stamp) < 0.6
    assert stamp.shape[0] == stamp.shape[1] == 1


def test_get_jitter_params():
    res = psf._get_jitter_params({})
    assert res["jitter_major"] > 0
    assert res["jitter_minor"] > 0
    assert np.isfinite(res["jitter_position_angle"])
    meta = {"jitter_major": 4}
    res = psf._get_jitter_params(meta)
    assert res["jitter_major"] == meta["jitter_major"]


def rms(stamp, coord):
    mu = np.sum(coord * stamp) / np.sum(stamp)
    sigma = np.sqrt(np.sum((coord - mu) ** 2 * stamp) / np.sum(stamp))
    return sigma


def test_evaluate_gaussian_fft():
    param = dict(jitter_major=8, jitter_minor=8, jitter_position_angle=0)
    shape = (19, 19)
    fft = psf._evaluate_gaussian_fft(param, shape, 0.008)
    stamp = np.fft.fftshift(np.fft.irfft2(fft, s=shape))
    tolerance = 0.01
    assert np.abs(np.sum(stamp) - 1) < tolerance

    # check that the gaussian's have the right rms's
    xx, yy = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    assert np.abs(rms(stamp, xx) - 1) < tolerance
    assert np.abs(rms(stamp, yy) - 1) < tolerance

    fft = psf._evaluate_gaussian_fft(param, shape, 0.004)
    stamp = np.fft.fftshift(np.fft.irfft2(fft, s=shape))
    assert np.abs(rms(stamp, xx) - 2) < tolerance
    assert np.abs(rms(stamp, yy) - 2) < tolerance

    # check that the position angle uses the right conventions
    param["jitter_major"] = 24
    param["jitter_minor"] = 8
    fft = psf._evaluate_gaussian_fft(param, shape, 0.008)
    stamp = np.fft.fftshift(np.fft.irfft2(fft, s=shape))
    assert np.abs(rms(stamp, yy) - 3) < tolerance
    assert np.abs(rms(stamp, xx) - 1) < tolerance

    # make sure we tilt to the right when the position angle
    # is mildly positive
    param["jitter_position_angle"] = 30
    fft = psf._evaluate_gaussian_fft(param, shape, 0.008)
    stamp = np.fft.fftshift(np.fft.irfft2(fft, s=shape))
    tophalf = np.s_[shape[0] // 2, :]
    assert np.sum((xx * stamp)[tophalf]) / np.sum(stamp[tophalf]) > 0


def test_add_jitter(setup_inputs):
    img = setup_inputs["image"]
    psf_ref_model = setup_inputs["psf_ref_model_f087"]

    # at least check that it runs and that the results are different?
    # it would be nice to check that the stamp sizes change in the
    # expected ways...
    img = deepcopy(img)
    psf_ref_model = deepcopy(psf_ref_model)
    psf_ref_model.meta.jitter_major = 8
    psf_ref_model.meta.jitter_minor = 8
    psf_ref_model.meta.jitter_position_angle = 0
    img.meta.guide_star.jitter_major = 16
    img.meta.guide_star.jitter_minor = 16
    img.meta.guide_star.jitter_position_angle = 0

    newstamps = psf.add_jitter(psf_ref_model, img)
    shape = psf_ref_model.psf.shape[-2:]
    npts = 5
    center = np.s_[
        shape[0] // 2 - npts : shape[0] // 2 + npts,
        shape[1] // 2 - npts : shape[1] // 2 + npts,
    ]
    idx = next(iter(np.ndindex(psf_ref_model.psf.shape[:-2])))
    xx, yy = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    assert rms(newstamps[idx][center], xx[center]) > rms(
        psf_ref_model.psf[idx][center], xx[center]
    )
    img.meta.guide_star.jitter_major = 0
    img.meta.guide_star.jitter_minor = 0
    newstamps = psf.add_jitter(psf_ref_model, img)
    assert rms(newstamps[idx][center], xx[center]) < rms(
        psf_ref_model.psf[idx][center], xx[center]
    )

# new functions
# _cart_to_polar
# _azimuthally_average_via_fft
# _downsample_by_interpolation
# create_convolution_kernel
# central_stamp

def test_azimuthally_average_via_fft(setup_inputs):
    psf_ref_model_f087 = setup_inputs["psf_ref_model_f087"]
    img = psf_ref_model_f087.psf[0, 0, 0].copy()
    img_avg = psf._azimuthally_average_via_fft(img, pixel_scale_ratio=0.5)
    img_cen = psf.central_stamp(img, size=3)
    img_avg_cen = psf.central_stamp(img_avg, size=3)
    oldstd = np.std(img_cen / img_cen[::-1, ::-1])
    newstd = np.std(img_avg_cen / img_avg_cen[::-1, ::-1])
    # in my tests oldstd is 0.09 and newstd is <1e-9.
    assert oldstd > 0.01
    assert newstd < 0.01

    grid_x, grid_y = np.mgrid[0:201, 0:201]
    gauss_model = Gaussian2D(1.0, 100, 100, 20, 20)
    gauss = gauss_model(grid_x, grid_y)
    smoothed = psf._azimuthally_average_via_fft(gauss)
    delta = gauss - smoothed

    assert np.mean(delta) < 1.0e-6
    assert mad_std(delta) < 1.0e-3




@pytest.mark.parametrize(
    "size",
    [40, 41],
)
def test_downsample_by_interpolation(size):
    pts = np.arange(size) - (size - 1) / 2
    xx, yy = np.meshgrid(pts, pts)
    gaussian = np.exp(-(xx / 5) ** 2 / 2 - (yy / 5) ** 2 / 2)

    gaussian_downsample = psf._downsample_by_interpolation(
        gaussian, downsample=4)
    assert gaussian_downsample.shape[0] < gaussian.shape[0] / 4 + 1

    def center(img):
        yy, xx = np.meshgrid(np.arange(img.shape[0]),
                             np.arange(img.shape[1]))
        xcen = np.sum(xx * img) / np.sum(img)
        ycen = np.sum(yy * img) / np.sum(img)
        return xcen, ycen

    cen = center(gaussian)
    cen_downsample = center(gaussian_downsample)
    cenfrac = cen[0] - gaussian.shape[1] / 2
    cen_downsamplefrac = cen_downsample[0] - gaussian_downsample.shape[1] / 2
    # test that both of these images are centered the same way
    assert np.abs(cenfrac - cen_downsamplefrac) < 0.0001


def test_create_convolution_kernel(setup_inputs):
    mod_f087 = setup_inputs["psf_ref_model_f087"]
    mod_f184 = setup_inputs["psf_ref_model_f184"]
    stamp_f087 = psf.central_stamp(mod_f087.psf[0, 1, 0], 91)
    stamp_f184 = psf.central_stamp(mod_f184.psf[0, 1, 0], 91)
    conv_kernel = psf.create_convolution_kernel(
        stamp_f087, stamp_f184)
    sz = 19
    diff_noconv = np.sum(psf.central_stamp(
        stamp_f087 - stamp_f184, sz) ** 2)
    mod_f184_conv = convolve(stamp_f087, conv_kernel)
    diff_conv = np.sum(psf.central_stamp(
        mod_f184_conv - stamp_f184, sz) ** 2)
    print(diff_conv, diff_noconv, diff_conv / diff_noconv)
    # FIXME
    # in my tests diff_conv is 1e-5, diff_noconv is 0.008, and the ratio
    # is 1.5e-3
    assert diff_conv < diff_noconv / 100


def test_central_stamp():
    img = np.random.randn(99, 99)
    cen = psf.central_stamp(img, 19)
    assert cen.shape[0] == 19
    cen = psf.central_stamp(img, 20)
    assert cen.shape[0] == 21  # needed to make it bigger to be central
    img = np.random.randn(40, 40)
    cen = psf.central_stamp(img, 20)
    assert cen.shape[0] == 20
    cen = psf.central_stamp(img, 21)
    assert cen.shape[0] == 22  # needed to make it bigger to be central
