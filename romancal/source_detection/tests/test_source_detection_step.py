"""
 Unit tests for the Roman source detection step code
"""

import os

import numpy as np
import pytest
from astropy import units as u
from astropy.convolution import Gaussian2DKernel
from roman_datamodels import maker_utils as testutil
from roman_datamodels.datamodels import ImageModel

from romancal.source_detection import SourceDetectionStep


@pytest.fixture
def setup_inputs():
    def _setup(nrows=100, ncols=100, noise=1.0):

        """Return ImageModel of lvl 2 image"""

        shape = (100, 100)  # size of test image
        wfi_image = testutil.mk_level2_image(shape=shape)
        wfi_image.data = u.Quantity(
            np.ones(shape, dtype=np.float32), u.electron / u.s, dtype=np.float32
        )
        wfi_image.meta.filename = "filename"

        # add noise to data
        if noise is not None:
            wfi_image.data += u.Quantity(
                noise * np.random.random(shape), u.electron / u.s, dtype=np.float32
            )

        # add dq array

        wfi_image.dq = np.zeros(shape, dtype=np.uint32)
        # construct ImageModel
        mod = ImageModel(wfi_image)

        return mod

    return _setup


# @pytest.mark.skipif(
#     os.environ.get("CI") == "true",
#     reason="Roman CRDS servers are not currently available outside the internal "
#     "network",
# )
def add_random_gauss(arr, x_positions, y_positions, min_amp=200, max_amp=500):

    """Add random 2D Gaussians to `arr` at specified positions,
    with random amplitudes from `min_amp` to  `max_amp`. Assumes
    units of e-/s."""

    # choosing a random seed for now, total randomness was causing issues
    np.random.seed(0)

    for i, x in enumerate(x_positions):
        y = y_positions[i]
        gauss = Gaussian2DKernel(2, x_size=21, y_size=21).array
        amp = np.random.randint(200, 700)
        arr[y - 10 : y + 11, x - 10 : x + 11] += (
            u.Quantity(gauss, u.electron / u.s, dtype=np.float32) * amp
        )


# @pytest.mark.skipif(
#     os.environ.get("CI") == "true",
#     reason="Roman CRDS servers are not currently available outside the internal "
#     "network",
# )
def test_source_detection_defaults(setup_inputs):

    """Test SourceDetectionStep with its default parameters. The detection
    threshold will be chosen based on the image's background level."""

    model = setup_inputs()

    # add in 12 sources, roughly evenly distributed
    # sort by true_x so they can be matched up to output

    true_x = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 82])
    true_y = np.array([26, 80, 44, 19, 66, 39, 29, 72, 54, 29, 80, 62])

    # at each position, place a 2d gaussian
    # randomly vary flux from 100 to 500 for each source
    add_random_gauss(model.data, true_x, true_y)

    # call SourceDetectionStep with default parameters
    sd = SourceDetectionStep()
    res = sd.process(model)

    # unpack output catalog array
    _, xcentroid, ycentroid, flux = res.meta.source_detection.tweakreg_catalog

    # sort based on x coordinate, like input
    ycentroid = [x for y, x in sorted(zip(xcentroid, ycentroid))]
    flux = [x for y, x in sorted(zip(xcentroid, flux))]
    xcentroid = sorted(xcentroid)

    # check that the number of input and output sources are the same
    assert len(flux) == len(true_x)

    # check that their locations agree
    # atol=0.2 seems to be the lowest safe value for this right now.

    assert np.allclose(np.abs(xcentroid - true_x), 0.0, atol=0.25)
    assert np.allclose(np.abs(ycentroid - true_y), 0.0, atol=0.25)


# @pytest.mark.skipif(
#     os.environ.get("CI") == "true",
#     reason="Roman CRDS servers are not currently available outside the internal "
#     "network",
# )
def test_source_detection_scalar_threshold(setup_inputs):

    """Test SourceDetectionStep using the option to choose a detection
    threshold for entire image."""

    model = setup_inputs()

    # add in 12 sources, roughly evenly distributed
    # sort by true_x so they can be matched up to output

    true_x = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 82])
    true_y = np.array([26, 80, 44, 19, 66, 39, 29, 72, 54, 29, 80, 62])

    # at each position, place a 2d gaussian
    # randomly vary flux from 100 to 500 for each source
    add_random_gauss(model.data, true_x, true_y)

    # call SourceDetectionStep with default parameters
    sd = SourceDetectionStep()
    sd.scalar_threshold = 2.0
    res = sd.process(model)

    # unpack output catalog array
    _, xcentroid, ycentroid, flux = res.meta.source_detection.tweakreg_catalog

    # sort based on x coordinate, like input
    ycentroid = [x for y, x in sorted(zip(xcentroid, ycentroid))]
    flux = [x for y, x in sorted(zip(xcentroid, flux))]
    xcentroid = sorted(xcentroid)

    # check that the number of input and output sources are the same
    assert len(flux) == len(true_x)

    # check that their locations agree
    # atol=0.2 seems to be the lowest safe value for this right now.

    assert np.allclose(np.abs(xcentroid - true_x), 0.0, atol=0.25)
    assert np.allclose(np.abs(ycentroid - true_y), 0.0, atol=0.25)


@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently available outside the internal "
    "network",
)
def test_outputs(tmp_path, setup_inputs):
    """Make sure `save_catalogs` and `output_cat_filetype` work correctly."""

    model = setup_inputs()

    # add a single source to image so a non-empty catalog is produced
    add_random_gauss(model.data, [50], [50])

    # run step and direct it to save catalog. default format should be asdf
    sd = SourceDetectionStep()
    sd.save_catalogs = True
    sd.process(model)
    # make sure file exists
    expected_output_path = os.path.join(tmp_path, "filename_tweakreg_catalog.asdf")
    assert os.path.isfile(expected_output_path)

    # run again, specifying that catalog should be saved in ecsv format
    sd = SourceDetectionStep()
    sd.save_catalogs = True
    sd.output_cat_filetype = "ecsv"
    expected_output_path = os.path.join(tmp_path, "filename_tweakreg_catalog.ecsv")
    assert os.path.isfile(expected_output_path)


# @pytest.mark.skipif(
#     os.environ.get("CI") == "true",
#     reason="Roman CRDS servers are not currently available outside the internal "
#     "network",
# )
def test_limiting_catalog_size(setup_inputs):

    """Test to make sure setting `max_sources` limits the size of the
    output catalog to contain only the N brightest sources"""

    model = setup_inputs()

    amps = [200, 300, 400]  # flux
    pos = [20, 50, 80]  # 3 sources in a line
    for i in range(3):
        xy = pos[i]
        gauss = Gaussian2DKernel(2, x_size=20, y_size=20).array
        model.data[xy - 10 : xy + 10, xy - 10 : xy + 10] += (
            u.Quantity(gauss, u.electron / u.s, dtype=np.float32) * amps[i]
        )

    sd = SourceDetectionStep()
    sd.max_sources = 2
    res = sd.process(model)

    _, xcentroid, ycentroid, flux = res.meta.source_detection.tweakreg_catalog

    # make sure only 2 of the three sources are returned in output catalog
    assert len(flux) == 2

    # and make sure one of them is not the first, dimmest source
    assert np.all(xcentroid > 22)  # dimmest position is really 20, give a
