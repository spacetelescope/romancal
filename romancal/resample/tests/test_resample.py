import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.table import QTable
from astropy.time import Time
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils
from roman_datamodels.maker_utils import mk_common_meta, mk_level2_image

from romancal.datamodels import ModelContainer
from romancal.lib.tests.helpers import word_precision_check
from romancal.resample import gwcs_drizzle, resample_utils
from romancal.resample.resample import (
    ResampleData,
    populate_mosaic_basic,
    populate_mosaic_individual,
)


# Helper function to create a mock input model with specified metadata
def create_mock_model(
    start_time,
    end_time,
    visit,
    segment,
    pass_,
    program,
    survey,
    optical_element,
    instrument_name,
):
    meta = mk_common_meta()
    mock_model = mk_level2_image(**{"meta": meta})
    mock_model.meta.exposure.start_time = Time(start_time, format="mjd")
    mock_model.meta.exposure.end_time = Time(end_time, format="mjd")
    mock_model.meta.exposure.mid_time = Time((start_time + end_time) / 2, format="mjd")
    mock_model.meta.observation.visit = visit
    mock_model.meta.observation.segment = segment
    mock_model.meta.observation["pass"] = pass_
    mock_model.meta.observation.program = program
    mock_model.meta.observation.survey = survey
    mock_model.meta.instrument.optical_element = optical_element
    mock_model.meta.instrument.name = instrument_name
    mock_model.meta.wcsinfo.vparity = -1
    mock_model.meta.wcsinfo.v3yangle = -60
    return mock_model


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
                    u.MJy / u.sr,
                    dtype=np.float32,
                ),
                "var_rnoise": u.Quantity(
                    rng.normal(1, 0.05, size=self.shape).astype(np.float32),
                    u.MJy**2 / u.sr**2,
                    dtype=np.float32,
                ),
                "var_poisson": u.Quantity(
                    rng.poisson(1, size=self.shape).astype(np.float32),
                    u.MJy**2 / u.sr**2,
                    dtype=np.float32,
                ),
                "var_flat": u.Quantity(
                    rng.uniform(0, 1, size=self.shape).astype(np.float32),
                    u.MJy**2 / u.sr**2,
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
    single = False
    blendheaders = False
    pixfrac = 0.8
    kernel = "turbo"
    fillval = 0.0
    wht_type = "exp"
    good_bits = "1"
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
    assert resample_data.good_bits == "0"
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
    assert output_models[0].meta.filename == resample_data.output_filename
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


@pytest.mark.parametrize("weight_type", ["ivm", "exptime"])
def test_resampledata_do_drizzle_default_single_exposure_weight_array(
    exposure_1,
    weight_type,
):
    """Test that resample methods return non-empty weight arrays."""

    input_models = ModelContainer(exposure_1)
    resample_data = ResampleData(input_models, wht_type=weight_type)

    output_models_many_to_one = resample_data.resample_many_to_one()
    output_models_many_to_many = resample_data.resample_many_to_many()

    assert np.any(output_models_many_to_one[0].weight > 0)
    assert np.any(output_models_many_to_many[0].weight > 0)


def test_populate_mosaic_basic_single_exposure(exposure_1):
    """
    Test the populate_mosaic_basic function with a given exposure.
    """
    input_models = ModelContainer(exposure_1)
    output_wcs = resample_utils.make_output_wcs(
        input_models,
        pscale_ratio=1,
        pscale=0.000031,
        rotation=0,
        shape=None,
        crpix=(0, 0),
        crval=(0, 0),
    )
    output_model = maker_utils.mk_datamodel(
        datamodels.MosaicModel, shape=tuple(output_wcs.array_shape)
    )

    populate_mosaic_basic(output_model, input_models=input_models)

    input_meta = [datamodel.meta for datamodel in input_models]

    assert output_model.meta.basic.time_first_mjd == np.min(
        [x.exposure.start_time.mjd for x in input_meta]
    )
    assert output_model.meta.basic.time_last_mjd == np.max(
        [x.exposure.end_time.mjd for x in input_meta]
    )
    assert output_model.meta.basic.time_mean_mjd == np.mean(
        [x.exposure.mid_time.mjd for x in input_meta]
    )
    assert output_model.meta.basic.visit == (
        input_meta[0].observation.visit
        if len({x.observation.visit for x in input_meta}) == 1
        else -1
    )
    assert output_model.meta.basic.segment == (
        input_meta[0].observation.segment
        if len({x.observation.segment for x in input_meta}) == 1
        else -1
    )
    assert output_model.meta.basic["pass"] == (
        input_meta[0].observation["pass"]
        if len({x.observation["pass"] for x in input_meta}) == 1
        else -1
    )
    assert output_model.meta.basic.program == (
        input_meta[0].observation.program
        if len({x.observation.program for x in input_meta}) == 1
        else "-1"
    )
    assert output_model.meta.basic.survey == (
        input_meta[0].observation.survey
        if len({x.observation.survey for x in input_meta}) == 1
        else "MULTIPLE"
    )
    assert (
        output_model.meta.basic.optical_element
        == input_meta[0].instrument.optical_element
    )
    assert output_model.meta.basic.instrument == input_meta[0].instrument.name


@pytest.mark.parametrize(
    "input_models_data, expected_output",
    [
        (
            [
                (
                    59000.0,
                    59000.5,
                    1,
                    1,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
                (
                    59000.5,
                    59001.0,
                    1,
                    1,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
            ],
            {
                "time_first_mjd": 59000.0,
                "time_last_mjd": 59001.0,
                "time_mean_mjd": 59000.5,
                "visit": 1,
                "segment": 1,
                "pass": 1,
                "program": "12345",
                "survey": "N/A",
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
        # different visits
        (
            [
                (
                    59000.0,
                    59000.5,
                    1,
                    1,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
                (
                    59000.5,
                    59001.0,
                    2,
                    1,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
            ],
            {
                "time_first_mjd": 59000.0,
                "time_last_mjd": 59001.0,
                "time_mean_mjd": 59000.5,
                "visit": -1,
                "segment": 1,
                "pass": 1,
                "program": "12345",
                "survey": "N/A",
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
        # different segments
        (
            [
                (
                    59000.0,
                    59000.5,
                    1,
                    1,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
                (
                    59000.5,
                    59001.0,
                    1,
                    2,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
            ],
            {
                "time_first_mjd": 59000.0,
                "time_last_mjd": 59001.0,
                "time_mean_mjd": 59000.5,
                "visit": 1,
                "segment": -1,
                "pass": 1,
                "program": "12345",
                "survey": "N/A",
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
        # different passes
        (
            [
                (
                    59000.0,
                    59000.5,
                    1,
                    1,
                    1,
                    "12345",
                    "HLS",
                    "F158",
                    "WFI",
                ),
                (
                    59000.5,
                    59001.0,
                    1,
                    1,
                    2,
                    "12345",
                    "EMS",
                    "F158",
                    "WFI",
                ),
            ],
            {
                "time_first_mjd": 59000.0,
                "time_last_mjd": 59001.0,
                "time_mean_mjd": 59000.5,
                "visit": 1,
                "segment": 1,
                "pass": -1,
                "program": "12345",
                "survey": "MULTIPLE",
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
        # different programs
        (
            [
                (
                    59000.0,
                    59000.5,
                    1,
                    1,
                    1,
                    "12345",
                    "N/A",
                    "F158",
                    "WFI",
                ),
                (
                    59000.5,
                    59001.0,
                    1,
                    1,
                    1,
                    "54321",
                    "N/A",
                    "F158",
                    "WFI",
                ),
            ],
            {
                "time_first_mjd": 59000.0,
                "time_last_mjd": 59001.0,
                "time_mean_mjd": 59000.5,
                "visit": 1,
                "segment": 1,
                "pass": 1,
                "program": "-1",
                "survey": "N/A",
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
        # different surveys
        (
            [
                (
                    59000.0,
                    59000.5,
                    1,
                    1,
                    1,
                    "12345",
                    "HLS",
                    "F158",
                    "WFI",
                ),
                (
                    59000.5,
                    59001.0,
                    1,
                    1,
                    1,
                    "12345",
                    "EMS",
                    "F158",
                    "WFI",
                ),
            ],
            {
                "time_first_mjd": 59000.0,
                "time_last_mjd": 59001.0,
                "time_mean_mjd": 59000.5,
                "visit": 1,
                "segment": 1,
                "pass": 1,
                "program": "12345",
                "survey": "MULTIPLE",
                "optical_element": "F158",
                "instrument": "WFI",
            },
        ),
    ],
)
def test_populate_mosaic_basic_different_observations(
    input_models_data, expected_output
):
    """Test that populate_mosaic_basic function works properly under different observational scenarios."""
    input_models = [create_mock_model(*data) for data in input_models_data]
    output_wcs = resample_utils.make_output_wcs(
        input_models,
        pscale_ratio=1,
        pscale=0.000031,
        rotation=0,
        shape=None,
        crpix=(0, 0),
        crval=(0, 0),
    )
    output_model = maker_utils.mk_datamodel(
        datamodels.MosaicModel, shape=tuple(output_wcs.array_shape)
    )

    # Act
    populate_mosaic_basic(output_model, input_models)

    # Assert
    assert output_model.meta.basic.time_first_mjd == expected_output["time_first_mjd"]
    assert output_model.meta.basic.time_last_mjd == expected_output["time_last_mjd"]
    assert output_model.meta.basic.time_mean_mjd == expected_output["time_mean_mjd"]
    assert output_model.meta.basic.visit == expected_output["visit"]
    assert output_model.meta.basic.segment == expected_output["segment"]
    assert output_model.meta.basic.program == expected_output["program"]
    assert output_model.meta.basic.survey == expected_output["survey"]
    assert output_model.meta.basic.optical_element == expected_output["optical_element"]
    assert output_model.meta.basic.instrument == expected_output["instrument"]


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
            "ra_corn1": 10.005109345783163,
            "dec_corn1": -0.001982743978690467,
            "ra_corn2": 10.006897960220385,
            "dec_corn2": 0.002676755917536623,
            "ra_corn3": 10.000733528663718,
            "dec_corn3": 0.005043059486913547,
            "ra_corn4": 9.998944914237953,
            "dec_corn4": 0.00038355958555111314,
            "orientat_local": 9.826983262839223,
            "orientat": 9.826978421513601,
            "projection": "TAN",
            "s_region": (
                "POLYGON ICRS 10.005109345783163 -0.001982743978690467 10.006897960220385 "
                "0.002676755917536623 10.000733528663718 0.005043059486913547 "
                "9.998944914237953 0.00038355958555111314 "
            ),
        }
    )

    input_models = ModelContainer(multiple_exposures)
    resample_data = ResampleData(input_models)

    output_model = resample_data.resample_many_to_one()[0]

    assert output_model.meta.wcsinfo.projection == expected.projection
    assert word_precision_check(output_model.meta.wcsinfo.s_region, expected.s_region)
    for key in expected.keys():
        if key not in ["projection", "s_region"]:
            assert np.allclose(output_model.meta.wcsinfo[key], expected[key])


def test_l3_individual_image_meta(multiple_exposures):
    """Test that the individual_image_meta is being populated"""
    input_models = ModelContainer(multiple_exposures)
    output_model = maker_utils.mk_datamodel(datamodels.MosaicModel)

    # Act
    populate_mosaic_individual(output_model, input_models)

    # Assert sizes are expected
    n_inputs = len(input_models)
    for value in output_model.meta.individual_image_meta.values():
        assert isinstance(value, QTable)
        assert len(value) == n_inputs

    # Assert spot check on filename, which is different for each mock input
    basic_table = output_model.meta.individual_image_meta.basic
    for idx, input in enumerate(input_models):
        assert input.meta.filename == basic_table["filename"][idx]
