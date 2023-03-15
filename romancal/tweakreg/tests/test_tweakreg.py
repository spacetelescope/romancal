from roman_datamodels import datamodels as rdm
from romancal.tweakreg.tweakreg_step import TweakRegStep, _common_name
from roman_datamodels import maker_utils
import os
import csv
import pytest
from astropy import units as u
from gwcs import coordinate_frames as cf
from astropy import coordinates as coord
from gwcs import wcs
import numpy as np
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from gwcs.geometry import SphericalToCartesian, CartesianToSpherical
import copy


def update_wcsinfo(input_dm):
    input_dm.meta.wcsinfo.v2_ref = 0.42955128282521254
    input_dm.meta.wcsinfo.v3_ref = -0.2479976768255853
    input_dm.meta.wcsinfo.vparity = -1
    input_dm.meta.wcsinfo.v3yangle = -999999
    input_dm.meta.wcsinfo.ra_ref = 270.0
    input_dm.meta.wcsinfo.dec_ref = 66.0
    input_dm.meta.wcsinfo.roll_ref = 60
    input_dm.meta.wcsinfo.s_region = (
        "POLYGON IRCS "
        "269.3318903230621 65.56866666048172 "
        "269.32578768154605 65.69246311613287 "
        "269.02457173246125 65.69201346248587 "
        "269.0333096074621 65.56870823657276 "
    )


def _create_tel2sky_model(input_dm):
    """Create the transform from telescope to sky.

    The transform is defined with a reference point in a Frame
    associated tih the telescope (V2, V3) in arcsec, the corresponding
    reference poiont on sky (RA_REF, DEC_REF) in deg, and the position angle
    at the center of the aperture, ROLL_REF in deg.
    """
    v2_ref = input_dm.meta.wcsinfo.v2_ref / 3600
    v3_ref = input_dm.meta.wcsinfo.v3_ref / 3600
    roll_ref = input_dm.meta.wcsinfo.roll_ref
    ra_ref = input_dm.meta.wcsinfo.ra_ref
    dec_ref = input_dm.meta.wcsinfo.dec_ref

    angles = np.array([v2_ref, -v3_ref, roll_ref, dec_ref, -ra_ref])
    axes = "zyxyz"
    rot = RotationSequence3D(angles, axes_order=axes)

    # The sky rotation expects values in deg.
    # This should be removed when models work with quantities.
    model = (
        (Scale(1 / 3600) & Scale(1 / 3600))
        | SphericalToCartesian(wrap_lon_at=180)
        | rot
        | CartesianToSpherical(wrap_lon_at=360)
    )
    model.name = "v23tosky"
    return model


def create_wcs(input_dm, shift_1=0, shift_2=0):
    # create a WCS object
    shape = input_dm.data.shape
    # plate scale
    # ps = 0.11

    # create necessary transformations
    distortion = Shift(-shift_1) & Shift(-shift_2)
    distortion.bounding_box = ((-0.5, shape[-1] - 0.5), (-0.5, shape[-2] - 0.5))
    tel2sky = _create_tel2sky_model(input_dm)

    # create required frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # create pipeline
    pipeline = [
        wcs.Step(detector, distortion),
        wcs.Step(v2v3, tel2sky),
        wcs.Step(world, None),
    ]

    wcs_obj = wcs.WCS(pipeline)

    input_dm.meta["wcs"] = wcs_obj


def create_base_image_source_catalog(tmpdir, output_filename):
    # Gaia sources within 55 arcsec of RA, Dec = (269.32817654010694, 65.57015584665392)
    # src_sky_coords = [
    #     (269.2552105901048, 65.59884065259847),
    #     (269.2651157350142, 65.59502199010211),
    #     (269.23673910887766, 65.59858116842837),
    #     (269.27657752347227, 65.60268442163424),
    #     (269.27126490063, 65.5841304128806),
    #     (269.260215971339, 65.58022019091797),
    #     (269.225691472433, 65.58523716834911),
    #     (269.2089808807485, 65.59633630661216),
    #     (269.27084535269614, 65.62119550470338),
    #     (269.2370651943156, 65.62049352703104),
    #     (269.2122283260543, 65.5844946603055),
    #     (269.2926624422311, 65.57930017167308),
    #     (269.2052826809533, 65.61143216607557),
    #     (269.223367873589, 65.5777887463616),
    #     (269.195799639922, 65.5986465884684),
    #     (269.2163600379664, 65.61854764297102),
    #     (269.19498681141954, 65.59305204030738),
    #     (269.3214620104345, 65.59761578448311),
    #     (269.205520837197, 65.61573891804704),
    #     (269.1976412371581, 65.61057091792554),
    #     (269.28545645442796, 65.62369217329329),
    #     (269.2046006670173, 65.61872095513272),
    #     (269.3292811705684, 65.59027646369461),
    #     (269.2305896971811, 65.56946754296247),
    # ]

    # create a temp CSV file to be used as source catalog
    header = ["x", "y"]
    src_detector_coords = [
        [1021.6395657386537, 1001.9680688036428],
        [889.6851249750171, 875.2777123675428],
        [1273.499439434086, 993.5038786380983],
        [727.8615629227143, 1129.1750570535874],
        [814.2934358141154, 513.8273842811739],
        [967.9577295602869, 383.9665173073404],
        [1434.9320548022406, 550.6737516239843],
        [1653.851347594757, 919.230502289874],
        [792.0611181707495, 1741.9254227109516],
        [1251.9046329516914, 1719.2685030386315],
        [1619.3448905312116, 526.029234069159],
        [526.4315475714666, 353.38542516962025],
        [1691.9161761028918, 1419.7621970386222],
        [1472.7883297482276, 303.13986296513644],
        [1831.8185677548852, 995.9783646194151],
        [1535.2674205660842, 1655.214971385401],
        [1847.615051489303, 810.303768966457],
        [121.24370005941364, 960.8651530803602],
        [1685.1763839680657, 1562.398437078787],
        [1796.8018161166435, 1391.3260046806265],
        [591.7697163030869, 1824.2134932247782],
        [1695.3052636025634, 1661.1395432452134],
        [20.1408570884073, 717.5507147948415],
        [1381.027514556914, 26.359513306689905],
    ]
    output = os.path.join(tmpdir, output_filename)
    with open(output, "w", encoding="UTF8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(src_detector_coords)


def add_tweakreg_catalog_attribute(tmpdir, input_dm):
    # create and add a mock source detection catalog
    tweakreg_catalog_filename = "base_image_sources.csv"
    create_base_image_source_catalog(tmpdir, tweakreg_catalog_filename)
    input_dm.meta["tweakreg_catalog"] = os.path.join(tmpdir, tweakreg_catalog_filename)
    return input_dm


@pytest.fixture
def base_image():
    def _base_image(shift_1=0, shift_2=0):
        l2 = maker_utils.mk_level2_image(shape=(2000, 2000))
        # update wcsinfo
        update_wcsinfo(l2)
        # add a dummy WCS object
        create_wcs(l2, shift_1=shift_1, shift_2=shift_2)
        l2_im = rdm.ImageModel(l2)
        return l2_im

    return _base_image


@pytest.mark.parametrize(
    "input, error_type",
    [
        (list(), ValueError),
        ([""], FileNotFoundError),
        ("", TypeError),
        ([1, 2, 3], TypeError),
    ],
)
def test_tweakreg_raises_error_on_invalid_input(input, error_type):
    """Test that TweakReg raises an error when an invalid input is provided."""
    with pytest.raises(Exception) as exec_info:
        TweakRegStep.call(input)

    assert type(exec_info.value) == error_type


def test_tweakreg_raises_attributeerror_on_missing_tweakreg_catalog(base_image):
    """Test that TweakReg raises an AttributeError if meta.tweakreg_catalog is missing."""
    img = base_image()
    with pytest.raises(Exception) as exec_info:
        TweakRegStep.call([img])

    assert type(exec_info.value) == AttributeError


def test_tweakreg_returns_modelcontainer(tmpdir, base_image):
    """Test that TweakReg returns a ModelContainer."""
    img = base_image()
    add_tweakreg_catalog_attribute(tmpdir, img)
    res = TweakRegStep.call([img])

    assert type(res) == rdm.ModelContainer


def test_tweakreg_updates_cal_step(tmpdir, base_image):
    """Test that TweakReg updates meta.cal_step with tweakreg = COMPLETE."""
    img = base_image()
    add_tweakreg_catalog_attribute(tmpdir, img)
    res = TweakRegStep.call([img])

    assert hasattr(res[0].meta.cal_step, "tweakreg")
    assert res[0].meta.cal_step.tweakreg == "COMPLETE"


def test_tweakreg_updates_group_id(tmpdir, base_image):
    """Test that TweakReg updates 'group_id' with a non-zero length string."""
    img = base_image()
    add_tweakreg_catalog_attribute(tmpdir, img)
    res = TweakRegStep.call([img])

    assert hasattr(res[0].meta, "group_id")
    assert len(res[0].meta.group_id) > 0


@pytest.mark.parametrize(
    "shift_1, shift_2, tolerance, is_small_correction",
    [
        (0, 0, 5, True),
        (0, 3, 5, True),
        (1, 1, 5, True),
        (5, 5, (5**2 + 5**2) ** 0.5, True),
        (5, 5, 5, False),
        (5, 5, 1, False),
        (5, 5, 3, False),
        (5, 10, 5, False),
    ],
)
def test_tweakreg_correction_magnitude(
    shift_1, shift_2, tolerance, is_small_correction, request
):
    """Test that TweakReg corrections are within tolerance. All the parametrized values are in arcsec."""
    img1 = request.getfixturevalue("base_image")()
    img2 = request.getfixturevalue("base_image")(shift_1=shift_1, shift_2=shift_2)
    img1_wcs = copy.deepcopy(img1.meta.wcs)
    img2_wcs = copy.deepcopy(img2.meta.wcs)

    step = TweakRegStep()
    step.tolerance = tolerance / 10.0

    assert step._is_wcs_correction_small(img1_wcs, img2_wcs) == is_small_correction


@pytest.mark.parametrize(
    "filename_list, expected_common_name",
    (
        (
            [
                "l1-sca1_cal.asdf",
                "l1-sca2_cal.asdf",
                "l1-sca3_cal.asdf",
            ],
            "l1-sca",
        ),
        (
            [
                "l1-270-66-gaia-2016-sca1_cal.asdf",
                "l1-270-66-gaia-2016-sca2_cal.asdf",
                "l1-270-66-gaia-2016-sca3_cal.asdf",
            ],
            "l1-270-66-gaia-2016-sca",
        ),
    ),
)
def test_tweakreg_common_name(filename_list, expected_common_name, request):
    """Test that TweakReg raises an error when an invalid input is provided."""
    img_list = []
    for filename in filename_list:
        img = request.getfixturevalue("base_image")()
        img.meta["filename"] = filename
        img_list.append(img)

    res = _common_name(img_list)

    assert res == expected_common_name


@pytest.mark.parametrize(
    "filename_list",
    (
        [
            "l1-sca1_cal.asdf",
            "l1-sca2_cal.asdf",
            "l1-sca3_cal.asdf",
        ],
        [
            "l1-270-66-gaia-2016-sca1_cal.asdf",
            "l1-270-66-gaia-2016-sca2_cal.asdf",
            "l1-270-66-gaia-2016-sca3_cal.asdf",
        ],
    ),
)
def test_tweakreg_common_name_raises_error_on_invalid_input(filename_list):
    """Test that TweakReg raises an error when an invalid input is provided."""
    with pytest.raises(Exception) as exec_info:
        _common_name(filename_list)

    assert type(exec_info.value) == TypeError


@pytest.mark.parametrize(
    "abs_refcat",
    (
        "GAIADR1",
        "GAIADR2",
        "GAIADR3",
    ),
)
def test_tweakreg_save_valid_abs_refcat(tmp_path, abs_refcat, request):
    """Test that TweakReg saves the catalog used for absolute astrometry."""
    img = request.getfixturevalue("base_image")()
    add_tweakreg_catalog_attribute(tmp_path, img)

    step = TweakRegStep()
    step.save_abs_catalog = True
    step.abs_refcat = abs_refcat

    step.process([img])

    # file will be written to this directory by default
    root_dir = request.config.rootdir

    assert os.path.exists(root_dir / f"fit_{abs_refcat.lower()}_ref.ecsv")
    # clean up
    os.remove(root_dir / f"fit_{abs_refcat.lower()}_ref.ecsv")


def test_tweakreg_raises_error_on_invalid_abs_refcat(tmp_path, base_image):
    """Test that TweakReg raises an error when an invalid abs_refcat is provided."""
    img = base_image()
    add_tweakreg_catalog_attribute(tmp_path, img)

    step = TweakRegStep()
    step.save_abs_catalog = True
    step.abs_refcat = "my_ref_cat"

    with pytest.raises(Exception) as exec_info:
        step.process([img])

    assert type(exec_info.value) == ValueError
