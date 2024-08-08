import os
from itertools import product

import astropy.units as u
import numpy as np
import pytest
from astropy import coordinates as coord
from astropy.modeling import models
from gwcs import coordinate_frames as cf
from gwcs import wcs as gwcs_wcs
from roman_datamodels.datamodels import ImageModel
from roman_datamodels.dqflags import pixel
from roman_datamodels.maker_utils import mk_level2_image, mk_sky_background

from romancal.datamodels import ModelLibrary
from romancal.skymatch import SkyMatchStep


def mk_gwcs(shape, sky_offset=[0, 0] * u.arcsec, rotate=0 * u.deg):
    # Example adapted from photutils:
    #   https://github.com/astropy/photutils/blob/
    #   2825356f1d876cacefb3a03d104a4c563065375f/photutils/datasets/make.py#L821
    rho = np.pi / 3.0
    # Roman plate scale:
    scale = (0.11 * u.arcsec / u.pixel).to_value(u.deg / u.pixel)

    shift_by_crpix = models.Shift((-shape[1] / 2) + 1) & models.Shift(
        (-shape[0] / 2) + 1
    )

    cd_matrix = np.array(
        [
            [-scale * np.cos(rho), scale * np.sin(rho)],
            [scale * np.sin(rho), scale * np.cos(rho)],
        ]
    )

    rotation = models.AffineTransformation2D(cd_matrix, translation=[0, 0])
    rotation.inverse = models.AffineTransformation2D(
        np.linalg.inv(cd_matrix), translation=[0, 0]
    )

    tan = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        197.8925 + sky_offset[0].to_value(u.deg),
        -1.36555556 + sky_offset[1].to_value(u.deg),
        180.0 + rotate.to_value(u.deg),
    )

    det2sky = shift_by_crpix | rotation | tan | celestial_rotation
    det2sky.name = "linear_transform"

    detector_frame = cf.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )

    sky_frame = cf.CelestialFrame(
        reference_frame=coord.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )

    pipeline = [(detector_frame, det2sky), (sky_frame, None)]

    return gwcs_wcs.WCS(pipeline)


def mk_image_model(
    rate_mean=0,
    rate_std=1e-4,
    sky_offset=[0, 0] * u.arcsec,
    rotation=0 * u.deg,
    image_shape=(100, 100),
    rng=np.random.default_rng(619),
):
    l2 = mk_level2_image(shape=image_shape)
    l2_im = ImageModel(l2)
    l2_im.data = u.Quantity(
        rng.normal(loc=rate_mean, scale=rate_std, size=l2_im.data.shape).astype(
            np.float32
        ),
        l2_im.data.unit,
    )

    l2_im.meta["wcs"] = mk_gwcs(image_shape, sky_offset=sky_offset, rotate=rotation)

    l2_im.meta["background"] = mk_sky_background(
        level=None, subtracted=False, method="None"
    )

    l2_im.meta.cal_step["skymatch"] = "INCOMPLETE"
    return l2_im


@pytest.fixture
def wfi_rate():
    return mk_image_model()


@pytest.fixture
def mk_sky_match_image_models():
    rng = np.random.default_rng(1)
    im1a = mk_image_model(rng=rng)
    im1b = mk_image_model(sky_offset=[2, 2] * u.arcsec)
    im2a = mk_image_model(rotation=30 * u.deg)
    im2b = mk_image_model(sky_offset=[4, 4] * u.arcsec, rotation=60 * u.deg)
    im3 = mk_image_model(rotation=60 * u.deg)

    # add "bad" data
    im1a, dq_mask = _add_bad_pixels(im1a, 1e6, 1e9)
    im1b, _ = _add_bad_pixels(im1b, 1e6, 1e9)
    im2a, _ = _add_bad_pixels(im2a, 5e6, 3e9)
    im2b, _ = _add_bad_pixels(im2b, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    return [im1a, im1b, im2a, im2b, im3], dq_mask


def _add_bad_pixels(im, sat_val, dont_use_val):
    # Add two types of "bad" pixels: 1) in the center of the image that will
    # lie in the intersection region common to all images (we assume images
    # are rotated around a common center) and 2) bad pixels at the corners
    # of the images that have a different flag and will be excluded from the
    # analysis of rotated images because they will lie outside of the
    # intersection region common to all images (and not because they are
    # masked out based on DQ array).
    mask = np.ones(im.data.shape, dtype=bool)
    # Add some "bad" pixels:
    # corners
    im_unit = im.data.unit
    im.data[:5, :5] = sat_val * im_unit
    im.data[-5:, :5] = sat_val * im_unit
    im.data[-5:, -5:] = sat_val * im_unit
    im.data[:5, -5:] = sat_val * im_unit

    im.dq[:5, :5] = pixel.SATURATED
    im.dq[-5:, :5] = pixel.SATURATED
    im.dq[-5:, -5:] = pixel.SATURATED
    im.dq[:5, -5:] = pixel.SATURATED

    mask[:5, :5] = False
    mask[-5:, :5] = False
    mask[-5:, -5:] = False
    mask[:5, -5:] = False

    cy, cx = (x // 2 for x in im.data.shape)
    cx -= 5
    cy -= 5

    # center
    im.data[cx : cx + 10, cy : cy + 10] = dont_use_val * im_unit
    im.dq[cx : cx + 10, cy : cy + 10] = pixel.DO_NOT_USE
    mask[cx : cx + 10, cy : cy + 10] = False

    return im, mask


@pytest.mark.parametrize(
    "skymethod, subtract, skystat, match_down",
    tuple(
        product(
            ["local", "match", "global", "global+match"],
            [False, True],
            ["median", "mean", "midpt", "mode"],
            [False, True],
        )
    ),
)
def test_skymatch(wfi_rate, skymethod, subtract, skystat, match_down):
    # test basic functionality and correctness of sky computations
    rng = np.random.default_rng(42)
    im1 = wfi_rate.copy()
    im2 = im1.copy()
    im3 = im1.copy()

    # add "bad" data
    im1, dq_mask = _add_bad_pixels(im1, 1e6, 1e9)
    im2, _ = _add_bad_pixels(im2, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    library = ModelLibrary([im1, im2, im3])

    # define some background:
    levels = [9.12, 8.28, 2.56]

    with library:
        for i, (im, lev) in enumerate(zip(library, levels)):
            im.data = rng.normal(loc=lev, scale=0.05, size=im.data.shape) * im.data.unit
            library.shelve(im, i)

    # exclude central DO_NOT_USE and corner SATURATED pixels
    result = SkyMatchStep.call(
        library,
        skymethod=skymethod,
        match_down=match_down,
        subtract=subtract,
        skystat=skystat,
        binwidth=0.2,
        nclip=0,
        dqbits="~DO_NOT_USE+SATURATED",
    )

    if skymethod in ["local", "global+match"]:
        ref_levels = levels

    elif skymethod == "match":
        lev0 = min(levels) if match_down else max(levels)
        ref_levels = np.array(levels) - lev0

    elif skymethod == "global":
        ref_levels = len(levels) * [min(levels)]

    sub_levels = np.array(levels) - np.array(ref_levels)

    with result:
        for i, (im, lev, rlev, slev) in enumerate(
            zip(result, levels, ref_levels, sub_levels)
        ):
            # check that meta was set correctly:
            assert im.meta.background.method == skymethod
            assert im.meta.background.subtracted == subtract

            # test computed/measured sky values if level is set:
            if not np.isclose(im.meta.background.level.value, 0):
                assert abs(im.meta.background.level.value - rlev) < 0.01

            # test
            if subtract:
                assert abs(np.mean(im.data[dq_mask]).value - slev) < 0.01
            else:
                assert abs(np.mean(im.data[dq_mask]).value - lev) < 0.01
            result.shelve(im, i, modify=False)


@pytest.mark.parametrize(
    "skymethod, subtract, skystat",
    tuple(product(["local", "match", "global"], [False, True], ["mean"])),
)
def test_skymatch_overlap(mk_sky_match_image_models, skymethod, subtract, skystat):
    # test that computations are performed only in the area of overlap
    # between images (bad pixels in the corners of rotated images are ignored).
    # Set 'nclip' to 0 in order to not clip bad pixels in computing mean.
    rng = np.random.default_rng(7)
    [im1a, im1b, im2a, im2b, im3], dq_mask = mk_sky_match_image_models

    library = ModelLibrary([im1a, im1b, im2a, im2b, im3])

    # define some background:
    levels = [9.12, 9.12, 8.28, 8.28, 2.56]

    with library:
        for i, (im, lev) in enumerate(zip(library, levels)):
            im.data = rng.normal(loc=lev, scale=0.01, size=im.data.shape) * im.data.unit
            library.shelve(im, i)

    # We do not exclude SATURATED pixels. They should be ignored because
    # images are rotated and SATURATED pixels in the corners are not in the
    # common intersection of all input images. This is the purpose of this test
    result = SkyMatchStep.call(
        library,
        skymethod=skymethod,
        match_down=True,
        subtract=subtract,
        skystat=skystat,
        nclip=0,
        dqbits="~DO_NOT_USE",  # specifically DO NOT add 'SATURATED' flag
    )

    if skymethod in ["local", "global+match"]:
        ref_levels = levels

    elif skymethod == "match":
        ref_levels = np.array(levels) - min(levels)

    elif skymethod == "global":
        ref_levels = len(levels) * [min(levels)]

    sub_levels = np.array(levels) - np.array(ref_levels)

    with result:
        for i, (im, lev, rlev, slev) in enumerate(
            zip(result, levels, ref_levels, sub_levels)
        ):
            # check that meta was set correctly:
            assert im.meta.background.method == skymethod
            assert im.meta.background.subtracted == subtract

            if skymethod in ["local", "global"]:
                # These two sky methods must fail because they do not take
                # into account (do not compute) overlap regions and use
                # entire images:
                assert abs(im.meta.background.level.value - rlev) < 0.1

                # test
                if subtract:
                    assert abs(np.mean(im.data[dq_mask]).value - slev) < 0.1
                else:
                    assert abs(np.mean(im.data[dq_mask]).value - lev) < 0.01
            else:
                # test computed/measured sky values if level is nonzero:
                if not np.isclose(im.meta.background.level.value, 0):
                    assert abs(im.meta.background.level.value - rlev) < 0.01

                # test
                if subtract:
                    assert abs(np.mean(im.data[dq_mask].value) - slev) < 0.01
                else:
                    assert abs(np.mean(im.data[dq_mask].value) - lev) < 0.01
            result.shelve(im, i, modify=False)


@pytest.mark.parametrize(
    "skymethod, subtract",
    tuple(product(["local", "match", "global", "global+match"], [False, True])),
)
def test_skymatch_2x(wfi_rate, skymethod, subtract):
    # Test that repetitive applications of skymatch produce the same results
    rng = np.random.default_rng(19)
    im1 = wfi_rate.copy()
    im2 = im1.copy()
    im3 = im1.copy()

    # add "bad" data
    im1, dq_mask = _add_bad_pixels(im1, 1e6, 1e9)
    im2, _ = _add_bad_pixels(im2, 5e6, 3e9)
    im3, _ = _add_bad_pixels(im3, 7e6, 1e8)

    library = ModelLibrary([im1, im2, im3])

    # define some background:
    levels = [9.12, 8.28, 2.56]

    with library:
        for i, (im, lev) in enumerate(zip(library, levels)):
            im.data = rng.normal(loc=lev, scale=0.05, size=im.data.shape) * im.data.unit
            library.shelve(im, i)

    # We do not exclude SATURATED pixels. They should be ignored because
    # images are rotated and SATURATED pixels in the corners are not in the
    # common intersection of all input images. This is the purpose of this test
    step = SkyMatchStep(
        skymethod=skymethod,
        match_down=True,
        subtract=subtract,
        skystat="mean",
        nclip=0,
        dqbits="~DO_NOT_USE+SATURATED",
    )

    result = step.run([im1, im2, im3])

    with result:
        model = result.borrow(0)
        assert model.meta.background.subtracted == step.subtract
        assert model.meta.background.level is not None
        result.shelve(model, 0, modify=False)

    # 2nd run.
    step.subtract = False
    result2 = step.run(result)

    with result2:
        model = result2.borrow(0)
        assert model.meta.background.subtracted == step.subtract
        assert model.meta.background.level is not None
        result2.shelve(model, 0, modify=False)

    # compute expected levels
    if skymethod in ["local", "global+match"]:
        ref_levels = levels

    elif skymethod == "match":
        ref_levels = np.array(levels) - min(levels)

    elif skymethod == "global":
        ref_levels = len(levels) * [min(levels)]

    sub_levels = np.array(levels) - np.array(ref_levels)

    # compare results
    with result2:
        for i, (im, lev, rlev, slev) in enumerate(
            zip(result2, levels, ref_levels, sub_levels)
        ):
            # check that meta was set correctly:
            assert im.meta.background.method == skymethod
            assert im.meta.background.subtracted == step.subtract

            # test computed/measured sky values:
            if not np.isclose(im.meta.background.level.value, 0, atol=1e-6):
                assert abs(im.meta.background.level.value - rlev) < 0.01

            # test
            if subtract:
                assert abs(np.mean(im.data[dq_mask]).value - slev) < 0.01
            else:
                assert abs(np.mean(im.data[dq_mask]).value - lev) < 0.01
            result2.shelve(im, i, modify=False)


@pytest.mark.parametrize(
    "input_type",
    [
        "ModelLibrary",
        "ASNFile",
        "DataModelList",
        "ASDFFilenameList",
    ],
)
def test_skymatch_always_returns_modellibrary_with_updated_datamodels(
    input_type,
    mk_sky_match_image_models,
    tmp_path,
    create_mock_asn_file,
):
    """Test that the SkyMatchStep always returns a ModelLibrary
    with updated data models after processing different input types."""

    os.chdir(tmp_path)
    [im1a, im1b, im2a, im2b, im3], dq_mask = mk_sky_match_image_models

    im1a.meta.filename = "im1a.asdf"
    im1b.meta.filename = "im1b.asdf"
    im2a.meta.filename = "im2a.asdf"
    im2b.meta.filename = "im2b.asdf"
    im3.meta.filename = "im3.asdf"

    library = ModelLibrary([im1a, im1b, im2a, im2b, im3])
    library._save(tmp_path)

    step_input_map = {
        "ModelLibrary": library,
        "ASNFile": create_mock_asn_file(
            tmp_path,
            members_mapping=[
                {"expname": im1a.meta.filename, "exptype": "science"},
                {"expname": im1b.meta.filename, "exptype": "science"},
                {"expname": im2a.meta.filename, "exptype": "science"},
                {"expname": im2b.meta.filename, "exptype": "science"},
                {"expname": im3.meta.filename, "exptype": "science"},
            ],
        ),
        "DataModelList": [im1a, im1b],
        "ASDFFilenameList": [im1a.meta.filename, im1b.meta.filename],
    }

    step_input = step_input_map.get(input_type)

    res = SkyMatchStep.call(step_input)

    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert model.meta.cal_step.skymatch == "COMPLETE"
            assert hasattr(model.meta, "background")
            res.shelve(model, i, modify=False)
