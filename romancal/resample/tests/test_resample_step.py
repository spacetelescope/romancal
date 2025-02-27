from functools import reduce

import numpy as np
import pytest
from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.table import QTable
from astropy.time import Time
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils

from romancal.datamodels import ModelLibrary
from romancal.resample import ResampleStep


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

    step = ResampleStep(good_bits=good_bits)
    res = step.run(img)

    assert res.meta.resample.good_bits == good_bits


def test_individual_image_meta(base_image):
    """Test that the individual_image_meta is being populated"""
    input_models = ModelLibrary([base_image() for _ in range(2)])
    output_model = ResampleStep().run(input_models)

    # Assert sizes are expected
    n_inputs = len(input_models)
    for value in output_model.meta.individual_image_meta.values():
        assert isinstance(value, QTable)
        assert len(value) == n_inputs

    # Assert spot check on filename, which is different for each mock input
    basic_table = output_model.meta.individual_image_meta.basic
    with input_models:
        for idx, input in enumerate(input_models):
            assert input.meta.filename == basic_table["filename"][idx]
            input_models.shelve(input, index=idx)

    assert "background" in output_model.meta.individual_image_meta


@pytest.mark.parametrize(
    "meta_overrides, expected_basic",
    [
        (  # 2 exposures, share visit, etc
            (
                {
                    "meta.observation.visit": 1,
                    "meta.observation.pass": 1,
                    "meta.observation.segment": 1,
                    "meta.observation.program": 1,
                    "meta.instrument.optical_element": "F158",
                    "meta.instrument.name": "WFI",
                },
                {
                    "meta.observation.visit": 1,
                    "meta.observation.pass": 1,
                    "meta.observation.segment": 1,
                    "meta.observation.program": 1,
                    "meta.instrument.optical_element": "F158",
                    "meta.instrument.name": "WFI",
                },
            ),
            {
                "visit": 1,
                "pass": 1,
                "segment": 1,
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
        (  # 2 exposures, different metadata
            (
                {
                    "meta.observation.visit": 1,
                    "meta.observation.pass": 1,
                    "meta.observation.segment": 1,
                    "meta.observation.program": 1,
                    "meta.instrument.optical_element": "F158",
                    "meta.instrument.name": "WFI",
                },
                {
                    "meta.observation.visit": 2,
                    "meta.observation.pass": 2,
                    "meta.observation.segment": 2,
                    "meta.observation.program": 2,
                    "meta.instrument.optical_element": "F062",
                    "meta.instrument.name": "WFI",
                },
            ),
            {
                "visit": 1,
                "pass": 1,
                "segment": 1,
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
    ],
)
def test_populate_mosaic_basic(base_image, meta_overrides, expected_basic):
    """Test that the basic mosaic metadata is being populated"""
    models = []
    for i, meta_override in enumerate(meta_overrides):
        model = base_image()

        model.meta.observation.observation_id = i

        model.meta.exposure.start_time = Time(59000 + i, format="mjd")
        model.meta.exposure.end_time = Time(59001 + i, format="mjd")
        model.meta.exposure.mid_time = Time(
            (model.meta.exposure.start_time.mjd + model.meta.exposure.end_time.mjd) / 2,
            format="mjd",
        )
        model.meta.exposure.exposure_time = 3600 * 24

        for key, value in meta_override.items():
            *parent_keys, child_key = key.split(".")
            setattr(reduce(getattr, parent_keys, model), child_key, value)
        models.append(model)

    input_models = ModelLibrary(models)
    output_model = ResampleStep().run(input_models)

    assert (
        output_model.meta.basic.time_first_mjd == models[0].meta.exposure.start_time.mjd
    )
    assert (
        output_model.meta.basic.time_last_mjd == models[-1].meta.exposure.end_time.mjd
    )
    assert output_model.meta.basic.time_mean_mjd == np.mean(
        [m.meta.exposure.mid_time.mjd for m in models]
    )

    for key, value in expected_basic.items():
        assert getattr(output_model.meta.basic, key) == value


@pytest.mark.parametrize(
    "input_pixel_area, pixel_scale_ratio, expected_pixel_area",
    [
        # (None, 1.0, None), # this cannot be tested since it causes the step to crash
        (1.0, 1.0, -999999.0),
        (1.0, 2.0, -999999.0 * 4.0),
    ],
)
def test_pixel_area_update(
    base_image, input_pixel_area, pixel_scale_ratio, expected_pixel_area
):
    # if input model has a non-None pixel resample should scale it by the square of the pixel_scale_ratio
    model = base_image()
    model.meta.photometry.pixel_area = input_pixel_area
    output_model = ResampleStep(pixel_scale_ratio=pixel_scale_ratio).run(
        ModelLibrary([model])
    )
    assert output_model.meta.photometry.pixel_area == expected_pixel_area
