import pytest
from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf
from roman_datamodels import datamodels

from romancal.assign_wcs.utils import add_s_region
from romancal.regtest import util
from romancal.resample import ResampleStep
from romancal.resample.l3_wcs import l3wcsinfo_to_wcs


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


def test_wcs_wcsinfo_matches(base_image):
    model = base_image()
    img = datamodels.ImageModel(model)
    add_s_region(img)

    output_model = ResampleStep().run(img)

    wcs_from_wcsinfo = l3wcsinfo_to_wcs(output_model.meta.wcsinfo)
    ra_mad, dec_mad = util.comp_wcs_grids_arcs(output_model.meta.wcs, wcs_from_wcsinfo)
    assert (ra_mad + dec_mad) / 2.0 < 1.0e-5


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
    model = base_image()
    model.meta.wcsinfo["vparity"] = -1
    model.meta.wcs.bounding_box = (-0.5, 99.5), (-0.5, 99.5)

    img = datamodels.ImageModel(model)
    add_s_region(img)

    img.data *= img.meta.photometry.conversion_megajanskys / img.data

    step = ResampleStep(good_bits=good_bits)
    res = step.run(img)

    assert res.meta.resample.good_bits == good_bits
