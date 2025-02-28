import asdf
import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.time import Time
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils

from romancal.assign_wcs.utils import add_s_region
from romancal.datamodels import ModelLibrary
from romancal.lib.tests.helpers import word_precision_check
from romancal.resample import ResampleStep, resample_utils


class WfiSca:
    def __init__(self, fiducial_world, pscale, shape, filename):
        self.fiducial_world = fiducial_world
        self.pscale = pscale
        self.shape = shape
        self.filename = filename

    def create_image(self):
        """
        Create a dummy L2 datamodel given the coordinates of the fiducial point,
        a pixel scale, and the image shape and filename.

        Returns
        -------
        datamodels.ImageModel
            An L2 ImageModel datamodel.
        """
        rng = np.random.default_rng(seed=13)
        l2 = maker_utils.mk_level2_image(
            shape=self.shape,
            **{
                "meta": {
                    "wcsinfo": {
                        "ra_ref": 10,
                        "dec_ref": 0,
                        "vparity": -1,
                        "v3yangle": -60,
                    },
                    "exposure": {
                        "exposure_time": 152.04000000000002,
                        "effective_exposure_time": 3.04 * 6 * 8,
                    },
                    "observation": {
                        "program": 5,
                        "execution_plan": 1,
                        "pass": 1,
                        "observation": 1,
                        "segment": 1,
                        "visit": 1,
                        "visit_file_group": 1,
                        "visit_file_sequence": 1,
                        "visit_file_activity": "01",
                        "exposure": 1,
                    },
                },
                "data": rng.poisson(2.5, size=self.shape).astype(np.float32),
                "var_rnoise": rng.normal(1, 0.05, size=self.shape).astype(np.float32),
                "var_poisson": rng.poisson(1, size=self.shape).astype(np.float32),
                "var_flat": rng.uniform(0, 1, size=self.shape).astype(np.float32),
            },
        )
        # data from WFISim simulation of SCA #01
        l2.meta.filename = self.filename
        l2.meta["wcs"] = create_wcs_object_without_distortion(
            fiducial_world=self.fiducial_world,
            pscale=self.pscale,
            shape=self.shape,
        )
        model = datamodels.ImageModel(l2)
        add_s_region(model)
        return model


def create_wcs_object_without_distortion(fiducial_world, pscale, shape):
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

    det2sky = shift | affine | scale | tan | celestial_rotation
    det2sky.name = "linear_transform"

    detector_frame = cf.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = cf.CelestialFrame(
        reference_frame=coord.FK5(), name="fk5", unit=(u.deg, u.deg)
    )

    pipeline = [(detector_frame, det2sky), (sky_frame, None)]

    wcs_obj = WCS(pipeline)

    wcs_obj.bounding_box = (
        (-0.5, shape[-1] - 0.5),
        (-0.5, shape[-2] - 0.5),
    )

    wcs_obj.pixel_shape = shape[::-1]
    wcs_obj.array_shape = shape

    return wcs_obj


@pytest.fixture
def wfi_sca1():
    sca = WfiSca(
        fiducial_world=(10, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0001_wfi01_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca2():
    sca = WfiSca(
        fiducial_world=(10.00139, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0001_wfi02_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca3():
    sca = WfiSca(
        fiducial_world=(10.00278, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0001_wfi03_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca4():
    sca = WfiSca(
        fiducial_world=(10, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0002_wfi01_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca5():
    sca = WfiSca(
        fiducial_world=(10.00139, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0002_wfi02_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca6():
    sca = WfiSca(
        fiducial_world=(10.00278, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0002_wfi03_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def exposure_1(wfi_sca1, wfi_sca2, wfi_sca3):
    """Returns a list with models corresponding to a dummy exposure 1."""
    # set the same exposure time for all SCAs
    for sca in [wfi_sca1, wfi_sca2, wfi_sca3]:
        sca.meta.exposure["start_time"] = Time(
            "2020-02-01T00:00:00", format="isot", scale="utc"
        )
        sca.meta.exposure["end_time"] = Time(
            "2020-02-01T00:02:30", format="isot", scale="utc"
        )
        sca.meta.observation.exposure = 1
        sca.meta.observation.observation_id = "1"
    return [wfi_sca1, wfi_sca2, wfi_sca3]


@pytest.fixture
def exposure_2(wfi_sca4, wfi_sca5, wfi_sca6):
    """Returns a list with models corresponding to a dummy exposure 2."""
    # set the same exposure time for all SCAs
    for sca in [wfi_sca4, wfi_sca5, wfi_sca6]:
        sca.meta.exposure["start_time"] = Time(
            "2020-05-01T00:00:00", format="isot", scale="utc"
        )
        sca.meta.exposure["end_time"] = Time(
            "2020-05-01T00:02:30", format="isot", scale="utc"
        )
        sca.meta.observation.exposure = 2
        sca.meta.observation.observation_id = "2"
    return [wfi_sca4, wfi_sca5, wfi_sca6]


@pytest.fixture
def multiple_exposures(exposure_1, exposure_2):
    """Returns a list with all the datamodels from exposure 1 and 2."""
    exposure_1.extend(exposure_2)
    return exposure_1


def get_resampled_wcs_pixel_scale(wcs):
    t = wcs.forward_transform
    for p in t.param_names:
        if p.startswith("factor"):
            return getattr(t, p).value


def test_resampledata_do_drizzle_many_to_one_default_no_rotation_single_exposure(
    exposure_1,
):
    """Test that output WCS encompass the entire combined input WCS region
    by checking that its extrema fall within the output WCS footprint.

    N.B.: since we are not providing the rotation parameter for the
    resample_utils.make_output_wcs method, the output WCS will have
    the same orientation (i.e. same PA) as the detector axes.
    """

    input_models = ModelLibrary(exposure_1)
    output_model = ResampleStep().run(input_models)

    output_min_value = np.min(output_model.meta.wcs.footprint())
    output_max_value = np.max(output_model.meta.wcs.footprint())

    def get_footprint(model, index):
        return model.meta.wcs.footprint()

    input_wcs_list = list(input_models.map_function(get_footprint, modify=False))

    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    # Assert
    np.testing.assert_array_less(output_min_value, expected_min_value)
    np.testing.assert_(output_max_value > expected_max_value)


def test_resampledata_do_drizzle_many_to_one_default_no_rotation_multiple_exposures(
    multiple_exposures,
):
    """Test that output WCS encompass the entire combined input WCS region
    by checking that its extrema fall within the output WCS footprint.

    N.B.: since we are not providing the rotation parameter for the
    resample_utils.make_output_wcs method, the output WCS will have
    the same orientation (i.e. same PA) as the detector axes.
    """

    input_models = ModelLibrary(multiple_exposures)
    output_model = ResampleStep().run(input_models)

    output_min_value = np.min(output_model.meta.wcs.footprint())
    output_max_value = np.max(output_model.meta.wcs.footprint())

    def get_footprint(model, index):
        return model.meta.wcs.footprint()

    input_wcs_list = list(input_models.map_function(get_footprint, modify=False))

    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    # Assert
    np.testing.assert_array_less(output_min_value, expected_min_value)
    np.testing.assert_(output_max_value > expected_max_value)


def test_resampledata_do_drizzle_many_to_one_default_rotation_0(exposure_1):
    """Test that output WCS encompass the entire combined input WCS region
    by checking that the output WCS footprint vertices are close to the
    expected vertices for the combined input WCS footprint.

    N.B.: in this case, rotation=0 will create a WCS that will be oriented North up.
    """

    input_models = ModelLibrary(exposure_1)

    output_model = ResampleStep(rotation=0).run(input_models)

    pscale = get_resampled_wcs_pixel_scale(output_model.meta.wcs)
    output_min_value = np.min(output_model.meta.wcs.footprint())
    output_max_value = np.max(output_model.meta.wcs.footprint())

    def get_footprint(model, index):
        return model.meta.wcs.footprint()

    input_wcs_list = list(input_models.map_function(get_footprint, modify=False))

    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    # Assert
    assert (
        (expected_min_value - 0.5001 * pscale) <= output_min_value <= expected_min_value
    )
    assert (
        (expected_max_value + 0.5001 * pscale) >= output_max_value >= expected_max_value
    )


def test_resampledata_do_drizzle_many_to_one_default_rotation_0_multiple_exposures(
    multiple_exposures,
):
    """Test that output WCS encompass the entire combined input WCS region
    by checking that the output WCS footprint vertices are close to the
    expected vertices for the combined input WCS footprint.

    N.B.: in this case, rotation=0 will create a WCS that will be oriented North up.
    """

    input_models = ModelLibrary(multiple_exposures)
    output_model = ResampleStep(rotation=0).run(input_models)

    output_min_value = np.min(output_model.meta.wcs.footprint())
    output_max_value = np.max(output_model.meta.wcs.footprint())
    pscale = get_resampled_wcs_pixel_scale(output_model.meta.wcs)

    def get_footprint(model, index):
        return model.meta.wcs.footprint()

    input_wcs_list = list(input_models.map_function(get_footprint, modify=False))

    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    # Assert
    assert (
        (expected_min_value - 0.5001 * pscale) <= output_min_value <= expected_min_value
    )
    assert (
        (expected_max_value + 0.5001 * pscale) >= output_max_value >= expected_max_value
    )


def test_resampledata_do_drizzle_many_to_one_single_input_model(wfi_sca1):
    """Test that the output of resample from a single input file creates a WCS
    footprint vertices that are close to the input WCS footprint's vertices."""

    input_models = ModelLibrary([wfi_sca1])
    output_model = ResampleStep(rotation=0).run(input_models)

    flat_1 = np.sort(wfi_sca1.meta.wcs.footprint().flatten())
    flat_2 = np.sort(output_model.meta.wcs.footprint().flatten())
    pscale = get_resampled_wcs_pixel_scale(output_model.meta.wcs)

    np.testing.assert_allclose(flat_1, flat_2, atol=0.5 * pscale)


def test_update_exposure_times_different_sca_same_exposure(exposure_1):
    """Test that update_exposure_times is properly updating the exposure parameters
    for a set of different SCAs belonging to the same exposure."""
    input_models = ModelLibrary(exposure_1)
    output_model = ResampleStep().run(input_models)

    # these three SCAs overlap, so the max exposure time is 3x.
    # get this time within 0.1 s.
    time_difference = (
        output_model.meta.resample.product_exposure_time
        - 3 * exposure_1[0].meta.exposure.effective_exposure_time
    )
    assert np.abs(time_difference) < 0.1
    assert (
        output_model.meta.basic.time_first_mjd
        == exposure_1[0].meta.exposure.start_time.mjd
    )
    assert (
        output_model.meta.basic.time_last_mjd
        == exposure_1[0].meta.exposure.end_time.mjd
    )


def test_update_exposure_times_same_sca_different_exposures(exposure_1, exposure_2):
    """Test that update_exposure_times is properly updating the exposure parameters
    for a set of the same SCA but belonging to different exposures."""
    input_models = ModelLibrary([exposure_1[0], exposure_2[0]])
    output_model = ResampleStep().run(input_models)

    with input_models:
        models = list(input_models)
        first_mjd = min(x.meta.exposure.start_time for x in models).mjd
        last_mjd = max(x.meta.exposure.end_time for x in models).mjd
        [input_models.shelve(model, i, modify=False) for i, model in enumerate(models)]

    # these exposures overlap perfectly so the max exposure time should
    # be equal to the individual time times two.
    time_difference = (
        output_model.meta.resample.product_exposure_time
        - 2 * exposure_1[0].meta.exposure.effective_exposure_time
    )
    assert np.abs(time_difference) < 0.1

    assert output_model.meta.basic.time_first_mjd == first_mjd

    assert output_model.meta.basic.time_last_mjd == last_mjd

    # likewise the per-pixel median exposure time is just 2x the individual
    # sca exposure time.
    time_difference = (
        output_model.meta.basic.max_exposure_time
        - 2 * exposure_1[0].meta.exposure.effective_exposure_time
    )
    assert np.abs(time_difference) < 0.1


def test_custom_wcs_input_small_overlap_no_rotation(wfi_sca1, wfi_sca3, tmp_path):
    """Test that resample can create a proper output in the edge case where the
    desired output WCS does not encompass the entire input datamodel but, instead, have
    just a small overlap."""
    input_models = ModelLibrary([wfi_sca1])
    wcs_path = tmp_path / "wcs.asdf"
    asdf.AsdfFile({"wcs": wfi_sca3.meta.wcs}).write_to(wcs_path)

    output_model = ResampleStep(output_wcs=str(wcs_path), rotation=0).run(input_models)

    np.testing.assert_allclose(output_model.meta.wcs(0, 0), wfi_sca3.meta.wcs(0, 0))


def test_custom_wcs_input_entire_field_no_rotation(multiple_exposures, tmp_path):
    """Test that resample can create a proper output that encompasses the entire
    combined FOV of the input datamodels."""
    input_models = ModelLibrary(multiple_exposures)

    with input_models:
        models = list(input_models)
        # create output WCS encompassing the entire exposure FOV
        output_wcs = resample_utils.make_output_wcs(
            models,
            rotation=0,
        )
        [input_models.shelve(model, i, modify=False) for i, model in enumerate(models)]

    wcs_path = tmp_path / "wcs.asdf"
    asdf.AsdfFile({"wcs": output_wcs}).write_to(wcs_path)
    output_model = ResampleStep(output_wcs=str(wcs_path)).run(input_models)

    output_min_value = np.min(output_model.meta.wcs.footprint())
    output_max_value = np.max(output_model.meta.wcs.footprint())
    pscale = get_resampled_wcs_pixel_scale(output_model.meta.wcs)

    def get_footprint(model, index):
        return model.meta.wcs.footprint()

    input_wcs_list = list(input_models.map_function(get_footprint, modify=False))

    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    assert (
        (expected_min_value - 0.5001 * pscale) <= output_min_value <= expected_min_value
    )
    assert (
        (expected_max_value + 0.5001 * pscale) >= output_max_value >= expected_max_value
    )


@pytest.mark.parametrize("weight_type", ["ivm", "exptime"])
def test_resampledata_do_drizzle_default_single_exposure_weight_array(
    exposure_1,
    weight_type,
):
    """Test that resample methods return non-empty weight arrays."""

    input_models = ModelLibrary(exposure_1)
    output_model = ResampleStep(weight_type=weight_type).run(input_models)
    assert np.any(output_model.weight > 0)


def test_l3_wcsinfo(multiple_exposures):
    """Test the population of the Level 3 wcsinfo block"""
    expected = maker_utils.mk_mosaic_wcsinfo(
        **{
            "ra_ref": 10.00292450000052,
            "dec_ref": 0.001534500000533253,
            "x_ref": 106.4579605214774,
            "y_ref": 80.66617532540977,
            "rotation_matrix": [
                [-0.9335804264969954, 0.3583679495458379],
                [0.3583679495458379, 0.9335804264969954],
            ],
            "pixel_scale": 3.100000000097307e-05,
            "pixel_scale_local": 3.099999999719185e-05,
            "pixel_shape": (161, 213),
            "ra_center": 10.002930353020417,
            "dec_center": 0.0015101325554100666,
            "ra_corn1": 10.005118261576513,
            "dec_corn1": -0.0020027691784169498,
            "ra_corn2": 10.006906876013732,
            "dec_corn2": 0.0026567307177480667,
            "ra_corn3": 10.000742444457124,
            "dec_corn3": 0.005023034287225611,
            "ra_corn4": 9.998953830031317,
            "dec_corn4": 0.00036353438578227396,
            "orientat_local": 20.999999978134802,
            "orientat": 20.99999999880985,
            "projection": "TAN",
            "s_region": (
                "POLYGON ICRS  10.005118262 -0.002002769 10.006906876 "
                "0.002656731 10.000742444 0.005023034 9.998953830 0.000363534"
            ),
        }
    )

    input_models = ModelLibrary(multiple_exposures)
    output_model = ResampleStep().run(input_models)

    assert output_model.meta.wcsinfo.projection == expected.projection
    assert word_precision_check(output_model.meta.wcsinfo.s_region, expected.s_region)
    for key in expected.keys():
        if key not in ["projection", "s_region"]:
            assert np.allclose(output_model.meta.wcsinfo[key], expected[key])
