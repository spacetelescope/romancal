import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.time import Time
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils

from romancal.datamodels import ModelContainer
from romancal.resample import gwcs_drizzle, resample_utils
from romancal.resample.resample import ResampleData


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
        rng = np.random.default_rng()
        l2 = maker_utils.mk_level2_image(
            shape=self.shape,
            **{
                "meta": {
                    "wcsinfo": {"ra_ref": 10, "dec_ref": 0, "vparity": -1},
                    "exposure": {
                        "exposure_time": 152.04000000000002,
                        "effective_exposure_time": 3.04 * 6 * 8,
                    },
                    "observation": {
                        "program": "00005",
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
                "data": u.Quantity(
                    rng.poisson(2.5, size=self.shape).astype(np.float32),
                    u.electron / u.s,
                    dtype=np.float32,
                ),
                "var_rnoise": u.Quantity(
                    rng.normal(1, 0.05, size=self.shape).astype(np.float32),
                    u.electron**2 / u.s**2,
                    dtype=np.float32,
                ),
                "var_poisson": u.Quantity(
                    rng.poisson(1, size=self.shape).astype(np.float32),
                    u.electron**2 / u.s**2,
                    dtype=np.float32,
                ),
                "var_flat": u.Quantity(
                    rng.uniform(0, 1, size=self.shape).astype(np.float32),
                    u.electron**2 / u.s**2,
                    dtype=np.float32,
                ),
            },
        )
        # data from WFISim simulation of SCA #01
        l2.meta.filename = self.filename
        l2.meta["wcs"] = create_wcs_object_without_distortion(
            fiducial_world=self.fiducial_world,
            pscale=self.pscale,
            shape=self.shape,
        )
        return datamodels.ImageModel(l2)


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
        filename="r0000501001001001001_01101_0001_WFI01_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca2():
    sca = WfiSca(
        fiducial_world=(10.00139, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0001_WFI02_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca3():
    sca = WfiSca(
        fiducial_world=(10.00278, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0001_WFI03_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca4():
    sca = WfiSca(
        fiducial_world=(10, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0002_WFI01_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca5():
    sca = WfiSca(
        fiducial_world=(10.00139, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0002_WFI02_cal.asdf",
    )

    return sca.create_image()


@pytest.fixture
def wfi_sca6():
    sca = WfiSca(
        fiducial_world=(10.00278, 0),
        pscale=(0.000031, 0.000031),
        shape=(100, 100),
        filename="r0000501001001001001_01101_0002_WFI03_cal.asdf",
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
        sca.meta.observation["exposure"] = 1
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
        sca.meta.observation["exposure"] = 2
    return [wfi_sca4, wfi_sca5, wfi_sca6]


@pytest.fixture
def multiple_exposures(exposure_1, exposure_2):
    """Returns a list with all the datamodels from exposure 1 and 2."""
    exposure_1.extend(exposure_2)
    return exposure_1


def test_resampledata_init(exposure_1):
    """Test that ResampleData can set initial values."""
    input_models = exposure_1
    output = "output.asdf"
    single = True
    blendheaders = False
    pixfrac = 0.8
    kernel = "turbo"
    fillval = 0.0
    wht_type = "exp"
    good_bits = 1
    pscale_ratio = 0.5
    pscale = 0.1
    kwargs = {"in_memory": False}

    resample_data = ResampleData(
        input_models,
        output=output,
        single=single,
        blendheaders=blendheaders,
        pixfrac=pixfrac,
        kernel=kernel,
        fillval=fillval,
        wht_type=wht_type,
        good_bits=good_bits,
        pscale_ratio=pscale_ratio,
        pscale=pscale,
        **kwargs,
    )

    # Assert
    assert resample_data.input_models == input_models
    assert resample_data.output_filename == output
    assert resample_data.pscale_ratio == pscale_ratio
    assert resample_data.single == single
    assert resample_data.blendheaders == blendheaders
    assert resample_data.pixfrac == pixfrac
    assert resample_data.kernel == kernel
    assert resample_data.fillval == fillval
    assert resample_data.weight_type == wht_type
    assert resample_data.good_bits == good_bits
    assert resample_data.in_memory == kwargs["in_memory"]


def test_resampledata_init_default(exposure_1):
    """Test instantiating ResampleData with default values."""
    input_models = exposure_1
    # Default parameter values

    resample_data = ResampleData(input_models)

    # Assert
    assert resample_data.input_models == input_models
    assert resample_data.output_filename is None
    assert resample_data.pscale_ratio == 1.0
    assert not resample_data.single
    assert resample_data.blendheaders
    assert resample_data.pixfrac == 1.0
    assert resample_data.kernel == "square"
    assert resample_data.fillval == "INDEF"
    assert resample_data.weight_type == "ivm"
    assert resample_data.good_bits == 0
    assert resample_data.in_memory


@pytest.mark.parametrize("input_models", [None, list(), [""], ModelContainer()])
def test_resampledata_init_invalid_input(input_models):
    """Test that ResampleData will raise an exception on invalid inputs."""
    with pytest.raises(Exception) as exec_info:
        ResampleData(input_models)

    assert type(exec_info.value) == ValueError


def test_resampledata_do_drizzle_many_to_one_default_no_rotation_single_exposure(
    exposure_1,
):
    """Test that output WCS encompass the entire combined input WCS region
    by checking that its extrema fall within the output WCS footprint.

    N.B.: since we are not providing the rotation parameter for the
    resample_utils.make_output_wcs method, the output WCS will have
    the same orientation (i.e. same PA) as the detector axes.
    """

    input_models = ModelContainer(exposure_1)
    resample_data = ResampleData(input_models)

    output_models = resample_data.resample_many_to_one()

    output_min_value = np.min(output_models[0].meta.wcs.footprint())
    output_max_value = np.max(output_models[0].meta.wcs.footprint())

    input_wcs_list = [sca.meta.wcs.footprint() for sca in input_models]
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

    input_models = ModelContainer(multiple_exposures)
    resample_data = ResampleData(input_models)

    output_models = resample_data.resample_many_to_one()

    output_min_value = np.min(output_models[0].meta.wcs.footprint())
    output_max_value = np.max(output_models[0].meta.wcs.footprint())

    input_wcs_list = [sca.meta.wcs.footprint() for sca in multiple_exposures]
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

    input_models = ModelContainer(exposure_1)
    resample_data = ResampleData(input_models, **{"rotation": 0})

    output_models = resample_data.resample_many_to_one()

    output_min_value = np.min(output_models[0].meta.wcs.footprint())
    output_max_value = np.max(output_models[0].meta.wcs.footprint())

    input_wcs_list = [sca.meta.wcs.footprint() for sca in exposure_1]
    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    # Assert
    np.testing.assert_allclose(output_min_value, expected_min_value)
    np.testing.assert_allclose(output_max_value, expected_max_value)


def test_resampledata_do_drizzle_many_to_one_default_rotation_0_multiple_exposures(
    multiple_exposures,
):
    """Test that output WCS encompass the entire combined input WCS region
    by checking that the output WCS footprint vertices are close to the
    expected vertices for the combined input WCS footprint.

    N.B.: in this case, rotation=0 will create a WCS that will be oriented North up.
    """

    input_models = ModelContainer(multiple_exposures)
    resample_data = ResampleData(input_models, **{"rotation": 0})

    output_models = resample_data.resample_many_to_one()

    output_min_value = np.min(output_models[0].meta.wcs.footprint())
    output_max_value = np.max(output_models[0].meta.wcs.footprint())

    input_wcs_list = [sca.meta.wcs.footprint() for sca in multiple_exposures]
    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    # Assert
    np.testing.assert_allclose(output_min_value, expected_min_value)
    np.testing.assert_allclose(output_max_value, expected_max_value)


def test_resampledata_do_drizzle_many_to_one_single_input_model(wfi_sca1):
    """Test that the output of resample from a single input file creates a WCS
    footprint vertices that are close to the input WCS footprint's vertices."""

    input_models = ModelContainer([wfi_sca1])
    resample_data = ResampleData(
        input_models, output=wfi_sca1.meta.filename, **{"rotation": 0}
    )

    output_models = resample_data.resample_many_to_one()

    flat_1 = np.sort(wfi_sca1.meta.wcs.footprint().flatten())
    flat_2 = np.sort(output_models[0].meta.wcs.footprint().flatten())

    # Assert
    assert len(output_models) == 1
    assert output_models[0].meta.basic.filename == resample_data.output_filename
    np.testing.assert_allclose(flat_1, flat_2)


def test_update_exposure_times_different_sca_same_exposure(exposure_1):
    """Test that update_exposure_times is properly updating the exposure parameters
    for a set of different SCAs belonging to the same exposure."""
    input_models = ModelContainer(exposure_1)
    resample_data = ResampleData(input_models)

    output_model = resample_data.resample_many_to_one()[0]

    exptime_tot = resample_data.resample_exposure_time(output_model)
    resample_data.update_exposure_times(output_model, exptime_tot)

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
    input_models = ModelContainer([exposure_1[0], exposure_2[0]])
    resample_data = ResampleData(input_models)

    output_model = resample_data.resample_many_to_one()[0]

    exptime_tot = resample_data.resample_exposure_time(output_model)
    resample_data.update_exposure_times(output_model, exptime_tot)

    assert len(resample_data.input_models.models_grouped) == 2

    # these exposures overlap perfectly so the max exposure time should
    # be equal to the individual time times two.
    time_difference = (
        output_model.meta.resample.product_exposure_time
        - 2 * exposure_1[0].meta.exposure.effective_exposure_time
    )
    assert np.abs(time_difference) < 0.1

    assert (
        output_model.meta.basic.time_first_mjd
        == min(x.meta.exposure.start_time for x in input_models).mjd
    )

    assert (
        output_model.meta.basic.time_last_mjd
        == max(x.meta.exposure.end_time for x in input_models).mjd
    )

    # likewise the per-pixel median exposure time is just 2x the individual
    # sca exposure time.
    time_difference = (
        output_model.meta.basic.max_exposure_time
        - 2 * exposure_1[0].meta.exposure.effective_exposure_time
    )
    assert np.abs(time_difference) < 0.1


@pytest.mark.parametrize(
    "name",
    ["var_rnoise", "var_poisson", "var_flat"],
)
def test_resample_variance_array(wfi_sca1, wfi_sca4, name):
    """Test that the mean value for the variance array lies within 1% of the
    expectation."""
    input_models = ModelContainer([wfi_sca1, wfi_sca4])
    resample_data = ResampleData(input_models, **{"rotation": 0})

    output_model = resample_data.blank_output.copy()
    output_model.meta["resample"] = {}
    driz = gwcs_drizzle.GWCSDrizzle(
        output_model,
        outwcs=resample_data.output_wcs,
        pixfrac=resample_data.pixfrac,
        kernel=resample_data.kernel,
        fillval=resample_data.fillval,
    )
    [driz.add_image(x.data, x.meta.wcs) for x in resample_data.input_models]

    resample_data.resample_variance_array(name, output_model)

    # combined variance is inversely proportional to the number of "measurements"
    expected_combined_variance_value = np.nanmean(
        [getattr(x, name) for x in input_models]
    ) / len(input_models)

    np.isclose(
        np.nanmean(getattr(output_model, name)).value,
        expected_combined_variance_value,
        atol=0.01,
    )


def test_custom_wcs_input_small_overlap_no_rotation(wfi_sca1, wfi_sca3):
    """Test that resample can create a proper output in the edge case where the
    desired output WCS does not encompass the entire input datamodel but, instead, have
    just a small overlap."""
    input_models = ModelContainer([wfi_sca1])
    resample_data = ResampleData(
        input_models,
        **{"output_wcs": wfi_sca3.meta.wcs, "rotation": 0},
    )

    output_models = resample_data.resample_many_to_one()

    np.testing.assert_allclose(output_models[0].meta.wcs(0, 0), wfi_sca3.meta.wcs(0, 0))


def test_custom_wcs_input_entire_field_no_rotation(multiple_exposures):
    """Test that resample can create a proper output that encompasses the entire
    combined FOV of the input datamodels."""
    input_models = ModelContainer(multiple_exposures)
    # create output WCS encompassing the entire exposure FOV
    output_wcs = resample_utils.make_output_wcs(
        input_models,
        rotation=0,
    )
    resample_data = ResampleData(
        input_models,
        **{"output_wcs": output_wcs},
    )

    output_models = resample_data.resample_many_to_one()

    output_min_value = np.min(output_models[0].meta.wcs.footprint())
    output_max_value = np.max(output_models[0].meta.wcs.footprint())

    input_wcs_list = [sca.meta.wcs.footprint() for sca in multiple_exposures]
    expected_min_value = np.min(np.stack(input_wcs_list))
    expected_max_value = np.max(np.stack(input_wcs_list))

    np.testing.assert_allclose(output_min_value, expected_min_value)
    np.testing.assert_allclose(output_max_value, expected_max_value)
