from io import BytesIO

import numpy as np
import pytest
from asdf import AsdfFile
from asdf_astropy.testing.helpers import assert_model_equal
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils

from romancal.datamodels import ModelLibrary
from romancal.resample import ResampleStep, resample_utils


class MockModel:
    def __init__(self, pixel_area, pixel_scale_ratio):
        self.meta = MockMeta(pixel_area, pixel_scale_ratio)


class MockMeta:
    def __init__(self, pixel_area, pixel_scale_ratio):
        self.photometry = MockPhotometry(pixel_area)
        self.resample = MockResample(pixel_scale_ratio)


class MockPhotometry:
    def __init__(self, pixel_area):
        self.pixel_area = pixel_area


class MockResample:
    def __init__(self, pixel_scale_ratio):
        self.pixel_scale_ratio = pixel_scale_ratio


class Mosaic:
    def __init__(self, fiducial_world, pscale, shape, filename, n_images):
        self.fiducial_world = fiducial_world
        self.pscale = pscale
        self.shape = shape
        self.filename = filename
        self.n_images = n_images

    def create_mosaic(self):
        """
        Create a dummy L3 datamodel given the coordinates of the fiducial point,
        a pixel scale, and the image shape and filename.

        Returns
        -------
        datamodels.MosaicModel
            An L3 MosaicModel datamodel.
        """
        l3 = maker_utils.mk_level3_mosaic(
            shape=self.shape,
            n_images=self.n_images,
        )
        # data from WFISim simulation of SCA #01
        l3.meta.filename = self.filename
        l3.meta["wcs"] = create_wcs_object_without_distortion(
            fiducial_world=self.fiducial_world,
            pscale=self.pscale,
            shape=self.shape,
        )
        # Call the forward transform so its value is pulled from disk
        _ = l3.meta.wcs.forward_transform
        return datamodels.MosaicModel(l3)


def create_wcs_object_without_distortion(fiducial_world, pscale, shape, **kwargs):
    """
    Create a simple WCS object without either distortion or rotation.

    Parameters
    ----------
    fiducial_world : tuple
        A pair of values corresponding to the fiducial's world coordinate.
    pscale : tuple
        A pair of values corresponding to the pixel scale in each axis.
    shape : tuple
        A pair of values specifying the dimensions of the WCS object.

    Returns
    -------
    gwcs.WCS
        A gwcs.WCS object.
    """
    # components of the model
    shift = models.Shift() & models.Shift()
    affine = models.AffineTransformation2D(
        matrix=[[1, 0], [0, 1]], translation=[0, 0], name="pc_rotation_matrix"
    )
    scale = models.Scale(pscale[0]) & models.Scale(pscale[1])
    tan = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        fiducial_world[0],
        fiducial_world[1],
        180,
    )

    # transforms between frames
    # detector -> sky
    det2sky = shift | affine | scale | tan | celestial_rotation
    det2sky.name = "linear_transform"

    # frames
    detector_frame = cf.Frame2D(
        name="detector",
        axes_order=(0, 1),
        axes_names=("x", "y"),
        unit=(u.pix, u.pix),
    )
    sky_frame = cf.CelestialFrame(
        reference_frame=coord.FK5(), name="fk5", unit=(u.deg, u.deg)
    )

    pipeline = [
        (detector_frame, det2sky),
        (sky_frame, None),
    ]

    wcs_obj = WCS(pipeline)

    wcs_obj.bounding_box = kwargs.get(
        "bounding_box",
        (
            (-0.5, shape[-1] - 0.5),
            (-0.5, shape[-2] - 0.5),
        ),
    )

    wcs_obj.pixel_shape = shape[::-1]
    wcs_obj.array_shape = shape

    return wcs_obj


@pytest.fixture
def asdf_wcs_file():
    def _create_asdf_wcs_file(tmp_path):
        file_path = tmp_path / "wcs.asdf"
        wcs_data = create_wcs_object_without_distortion(
            (10, 0),
            (0.000031, 0.000031),
            (100, 100),
        )
        wcs = {"wcs": wcs_data}
        with AsdfFile(wcs) as af:
            af.write_to(file_path)
        return str(file_path)

    return _create_asdf_wcs_file


@pytest.mark.parametrize(
    "vals, name, min_vals, expected",
    [
        ([1, 2], "list1", None, [1, 2]),
        ([None, None], "list2", None, None),
        ([1, 2], "list4", [0, 0], [1, 2]),
    ],
)
def test_check_list_pars_valid(vals, name, min_vals, expected):
    step = ResampleStep()

    result = step._check_list_pars(vals, name, min_vals)
    assert result == expected


def test_load_custom_wcs_no_file():
    step = ResampleStep()
    result = step._load_custom_wcs(None, (512, 512))
    assert result is None


def test_load_custom_wcs_missing_output_shape(asdf_wcs_file):
    with pytest.raises(ValueError):
        step = ResampleStep()
        step._load_custom_wcs(asdf_wcs_file, None)


def test_load_custom_wcs_invalid_file(tmp_path):
    step = ResampleStep()
    invalid_file = tmp_path / "invalid.asdf"
    with open(invalid_file, "w") as f:
        f.write("invalid asdf file")

    with pytest.raises(ValueError):
        step._load_custom_wcs(str(invalid_file), (512, 512))


def test_load_custom_wcs_asdf_without_wcs_attribute(tmp_path):
    step = ResampleStep()
    file_path = tmp_path / "asdf_file.asdf"
    wcs = {}
    with AsdfFile(wcs) as af:
        af.write_to(file_path)

    with pytest.raises(KeyError):
        step._load_custom_wcs(str(file_path), (100, 100))


def _assert_frame_equal(a, b):
    """ Copied from `gwcs`'s test_wcs.py """
    __tracebackhide__ = True

    assert type(a) is type(b)

    if a is None:
        return

    if not isinstance(a, cf.CoordinateFrame):
        return a == b

    assert a.name == b.name  # nosec
    assert a.axes_order == b.axes_order  # nosec
    assert a.axes_names == b.axes_names  # nosec
    assert a.unit == b.unit  # nosec
    assert a.reference_frame == b.reference_frame  # nosec


def _assert_wcs_equal(a, b):
    """ Based on corresponding function from `gwcs`'s test_wcs.py """
    assert a.name == b.name  # nosec

    assert a.pixel_shape == b.pixel_shape
    assert a.array_shape == b.array_shape
    if a.array_shape is not None:
        assert a.array_shape == b.pixel_shape[::-1]

    assert len(a.available_frames) == len(b.available_frames)  # nosec
    for a_step, b_step in zip(a.pipeline, b.pipeline, strict=False):
        _assert_frame_equal(a_step.frame, b_step.frame)
        assert_model_equal(a_step.transform, b_step.transform)


@pytest.mark.parametrize(
    "wcs_type", ["obj", "bytes"],
)
def test_load_custom_wcs_obj(exposure_1, wcs_type):
    shape = (123, 456)
    wcs_obj = create_wcs_object_without_distortion(
        (10, 0),
        (0.000031, 0.000031),
        shape,
    )
    step = ResampleStep()

    if wcs_type == "obj":
        step.output_wcs = wcs_obj

    elif wcs_type == "bytes":
        # write to a byte stream:
        bstream = BytesIO()
        AsdfFile({"wcs": wcs_obj}).write_to(bstream)
        bstream.seek(0)
        step.output_wcs = bstream

    result = step.process(ModelLibrary(exposure_1))
    if wcs_type == "bytes":
        bstream.close()

    assert result["data"].shape == shape

    outwcs = result["meta"]["wcs"]
    assert outwcs is not wcs_obj
    _assert_wcs_equal(outwcs, wcs_obj)


@pytest.mark.parametrize(
    "vals, name, min_vals, expected",
    [
        ([1, 2], "test", [0, 0], [1, 2]),
        ([None, None], "test", [0, 0], None),
        ([0, 0], "test", [0, 0], [0, 0]),
        ([1, 1], "test", [0, 0], [1, 1]),
        ([0, 1], "test", [0, 0], [0, 1]),
        ([1, 0], "test", [0, 0], [1, 0]),
    ],
)
def test_check_list_pars(vals, name, min_vals, expected):
    step = ResampleStep()

    result = step._check_list_pars(vals, name, min_vals)
    assert result == expected


@pytest.mark.parametrize(
    "vals, name, min_vals",
    [
        ([None, 2], "test", [0, 0]),
        ([1, None], "test", [0, 0]),
        ([1], "test", [0, 0]),
        ([1, 2, 3], "test", [0, 0]),
        ([None, None, None], "test", [0, 0]),
        ([1, 2], "test", [2, 2]),
    ],
)
def test_check_list_pars_exception(vals, name, min_vals):
    step = ResampleStep()

    with pytest.raises(ValueError):
        step._check_list_pars(vals, name, min_vals)


@pytest.mark.parametrize(
    """pixel_area, pixel_scale_ratio,
    expected_steradians, expected_arcsecsq""",
    [
        # Happy path tests
        (1.0, 2.0, 4.0, 4.0),
        (2.0, 0.5, 0.5, 0.5),
        (0.0, 2.0, 0.0, 0.0),
        (1.0, 0.0, 0.0, 0.0),
        (None, 2.0, None, 4.0),
        (1.0, 2.0, 4.0, None),
    ],
)
def test_update_phot_keywords(
    pixel_area,
    pixel_scale_ratio,
    expected_steradians,
    expected_arcsecsq,
):
    step = ResampleStep()
    model = MockModel(pixel_area, pixel_scale_ratio)

    step.update_phot_keywords(model)

    assert model.meta.photometry.pixel_area == expected_steradians


@pytest.mark.parametrize(
    "good_bits, dq_array, expected_output",
    [
        (
            "~DO_NOT_USE+NON_SCIENCE",
            np.array([[0, 2**0, 0], [0, 2**13, 0], [0, 2**9, 0]]),
            np.array([[1, 0, 1], [1, 1, 1], [1, 0, 1]]),
        ),
        (
            "~513",
            np.array([[0, 2**0, 0], [2**13, 0, 0], [0, 0, 2**9]]),
            np.array([[1, 0, 1], [1, 1, 1], [1, 1, 0]]),
        ),
        (
            "~1+512",
            np.array([[2**13, 0, 0], [0, 0, 2**9], [0, 2**0, 0]]),
            np.array([[1, 1, 1], [1, 1, 0], [1, 0, 1]]),
        ),
        (
            "~1,512",
            np.array([[2**13, 2**0, 2**9], [0, 0, 0], [0, 0, 0]]),
            np.array([[1, 0, 0], [1, 1, 1], [1, 1, 1]]),
        ),
        (
            "LOW_QE+NONLINEAR",
            np.array([[2**13, 2**0, 2**16], [0, 0, 0], [0, 0, 0]]),
            np.array([[1, 0, 1], [1, 1, 1], [1, 1, 1]]),
        ),
        (
            "73728",
            np.array([[0, 0, 0], [0, 2**13, 2**0], [0, 2**16, 0]]),
            np.array([[1, 1, 1], [1, 1, 0], [1, 1, 1]]),
        ),
        (
            "8192+65536",
            np.array([[0, 0, 0], [0, 2**13, 0], [2**0, 2**16, 0]]),
            np.array([[1, 1, 1], [1, 1, 1], [0, 1, 1]]),
        ),
        (
            "8192,65536",
            np.array([[0, 0, 0], [0, 2**13, 0], [0, 2**16, 2**0]]),
            np.array([[1, 1, 1], [1, 1, 1], [1, 1, 0]]),
        ),
    ],
)
def test_build_driz_weight_multiple_good_bits(
    base_image, good_bits, dq_array, expected_output
):
    data_shape = dq_array.shape
    img1 = base_image()
    img1.dq = dq_array
    img1.data = np.ones(data_shape)

    result = resample_utils.build_driz_weight(
        img1, weight_type=None, good_bits=good_bits
    )

    np.testing.assert_array_equal(result, expected_output)


@pytest.mark.parametrize(
    "good_bits",
    [
        "~DO_NOT_USE+NON_SCIENCE",
        "~513",
        "~1+512",
        "~1,512",
        "LOW_QE+NONLINEAR",
        "73728",
        "8192+65536",
        "8192,65536",
    ],
)
def test_set_good_bits_in_resample_meta(base_image, good_bits):
    model = maker_utils.mk_level2_image(shape=(100, 100))
    model.meta.wcsinfo["vparity"] = -1

    img = datamodels.ImageModel(model)

    img.data *= img.meta.photometry.conversion_megajanskys / img.data

    step = ResampleStep

    res = step.call(img, good_bits=good_bits)

    assert res.meta.resample.good_bits == good_bits


@pytest.mark.parametrize("weight_type", ["ivm", "exptime", None])
def test_build_driz_weight_different_weight_type(base_image, weight_type):
    rng = np.random.default_rng()
    img1 = base_image()
    # update attributes that will be used in building the weight array
    img1.meta.exposure.exposure_time = 10
    img1.var_rnoise = rng.normal(1, 0.1, size=img1.shape)
    # build the drizzle weight array
    result = resample_utils.build_driz_weight(
        img1, weight_type=weight_type, good_bits="~DO_NOT_USE+NON_SCIENCE"
    )

    expected_results = {
        "ivm": img1.var_rnoise**-1,
        "exptime": np.ones(img1.shape, dtype=img1.data.dtype)
        * img1.meta.exposure.exposure_time,
        None: np.ones(img1.shape, dtype=img1.data.dtype),
    }

    np.testing.assert_array_almost_equal(expected_results.get(weight_type), result)
