import pytest
from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils

from romancal.resample import ResampleStep


class MockModel:
    def __init__(self, pixelarea_steradians, pixelarea_arcsecsq, pixel_scale_ratio):
        self.meta = MockMeta(
            pixelarea_steradians, pixelarea_arcsecsq, pixel_scale_ratio
        )


class MockMeta:
    def __init__(self, pixelarea_steradians, pixelarea_arcsecsq, pixel_scale_ratio):
        self.photometry = MockPhotometry(pixelarea_steradians, pixelarea_arcsecsq)
        self.resample = MockResample(pixel_scale_ratio)


class MockPhotometry:
    def __init__(self, pixelarea_steradians, pixelarea_arcsecsq):
        self.pixelarea_steradians = pixelarea_steradians
        self.pixelarea_arcsecsq = pixelarea_arcsecsq


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
        l3.meta.wcs.forward_transform
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
    """pixelarea_steradians, pixelarea_arcsecsq, pixel_scale_ratio,
    expected_steradians, expected_arcsecsq""",
    [
        # Happy path tests
        (1.0, 1.0, 2.0, 4.0, 4.0),
        (2.0, 2.0, 0.5, 0.5, 0.5),
        (0.0, 0.0, 2.0, 0.0, 0.0),
        (1.0, 1.0, 0.0, 0.0, 0.0),
        (None, 1.0, 2.0, None, 4.0),
        (1.0, None, 2.0, 4.0, None),
    ],
)
def test_update_phot_keywords(
    pixelarea_steradians,
    pixelarea_arcsecsq,
    pixel_scale_ratio,
    expected_steradians,
    expected_arcsecsq,
):
    step = ResampleStep()
    model = MockModel(pixelarea_steradians, pixelarea_arcsecsq, pixel_scale_ratio)

    step.update_phot_keywords(model)

    assert model.meta.photometry.pixelarea_steradians == expected_steradians
    assert model.meta.photometry.pixelarea_arcsecsq == expected_arcsecsq
