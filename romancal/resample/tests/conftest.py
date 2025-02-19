import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.time import Time
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels, maker_utils


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
