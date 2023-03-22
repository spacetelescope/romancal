import copy
import csv
import os

import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from numpy.testing import assert_allclose
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.tweakreg.astrometric_utils import compute_radius, get_catalog
from romancal.tweakreg.tweakreg_step import TweakRegStep, _common_name


def update_wcsinfo(input_dm):
    """
    Update WCSInfo with realistic data (i.e. obtained from romanisim simulations).

    Parameters
    ----------
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.
    """
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
    """
    Create the transform from telescope to sky.

    The transform is defined with a reference point in a Frame
    associated with the telescope (V2, V3) in arcsec, the corresponding
    reference point on sky (RA_REF, DEC_REF) in deg, and the position angle
    at the center of the aperture, ROLL_REF in deg.

    Parameters
    ----------
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.

    Returns
    -------
    astropy.modeling.core.CompoundModel
        An astropy compound model with the transformations to map
        the telescope reference system onto the world reference systems.
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
        (Scale(0.1 / 3600) & Scale(0.1 / 3600))
        | SphericalToCartesian(wrap_lon_at=180)
        | rot
        | CartesianToSpherical(wrap_lon_at=360)
    )
    model.name = "v23tosky"
    return model


def create_basic_wcs(
    img_shape: tuple = (100, 100),
    ref_pix: tuple = (0, 0),
    ref_val: tuple[u.Quantity, u.Quantity] = (
        u.Quantity("10 deg"),
        u.Quantity("0 deg"),
    ),
    pix_scale: tuple[u.Quantity, u.Quantity] = (
        u.Quantity("0.1 arcsec"),
        u.Quantity("0.1 arcsec"),
    ),
    theta: u.Quantity = u.Quantity("0 deg"),
):
    """
    Creates a basic WCS (no distortion) to map pixel coordinates
    onto celestial spherical coordinates.

    Parameters
    ----------
    img_shape : tuple, optional
        Dimension of the pixel coordinates, by default (100, 100).
    ref_pix : tuple, optional
        Projection plane coordinates of the reference pixel, by default (0, 0)
    ref_val : tuple[u.Quantity, u.Quantity], optional
        Celestial coordinates at the reference pixel in degrees,
        by default (10, 0)*u.degree
    pix_scale : tuple[u.Quantity, u.Quantity], optional
        Size of a pixel side in arcsec, by default (0.1, 0.1)*u.arcsec
    theta : u.Quantity, optional
        CCW rotation angle of the projection plane, by default 0*u.deg

    Returns
    -------
    type: gwcs.wcs.WCS object
        A gwcs.wcs.WCS object with a model f : (x, y) |-> (ra, dec).

    Notes
    -----
    With the default parameters, the result will be a WCS for an image
    that's 100x100 pixels (@0.1"/pix in both axis) and reference
    coordinates (ra=10 deg, dec=0 deg, ICRS) sitting at the origin of the
    pixel coordinates frame (0,0).

    The WCS object created by this method generates a pipeline that
    does not contain the required steps to validate against the
    TweakReg pipeline.
    """

    # linear transformations
    shift_pixel_coords = models.Shift(-ref_pix[0]) & models.Shift(-ref_pix[1])
    scaling_matrix = np.array(
        [[pix_scale.to("deg")[0].value, 0], [0, pix_scale.to("deg")[1].value]]
    )
    rot_matrix = np.array(
        [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]
    )
    matrix = np.dot(rot_matrix, scaling_matrix)
    affine_transf = models.AffineTransformation2D(matrix)
    affine_transf.inverse = models.AffineTransformation2D(np.linalg.inv(matrix))
    tan = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        ref_val[0].value, ref_val[1].value, 180
    )
    det2sky = shift_pixel_coords | affine_transf | tan | celestial_rotation
    det2sky.name = "linear_transform"

    # frames
    detector_frame = cf.Frame2D(
        name="detector", unit=(u.pix, u.pix), axes_names=("x", "y")
    )
    cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    sky_frame = cf.CelestialFrame(
        name="sky", unit=(u.deg, u.deg), reference_frame=coord.ICRS()
    )

    pipeline = [(detector_frame, det2sky), (sky_frame, None)]
    wcsobj = wcs.WCS(pipeline)
    wcsobj.bounding_box = ((-0.5, img_shape[0] + 0.5), (-0.5, img_shape[1] + 0.5))

    return wcsobj


def create_wcs_for_tweakreg_pipeline(input_dm, shift_1=0, shift_2=0):
    """
    Create a basic WCS object (with optional shift) and
    append it to the input_dm.meta attribute.

    The WCS object will have a pipeline with the required
    steps to validate against the TweakReg pipeline.

    Parameters
    ----------
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.
    shift_1 : int, optional
        Shift to be applied in the direction of the first axis, by default 0
    shift_2 : int, optional
        Shift to be applied in the direction of the second axis, by default 0
    """

    shape = input_dm.data.shape

    # create necessary transformations
    distortion = Shift(-shift_1) & Shift(-shift_2)
    distortion.bounding_box = ((-0.5, shape[-1] + 0.5), (-0.5, shape[-2] + 0.5))
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


def create_base_image_source_catalog(tmp_path, output_filename, catalog_data=None):
    """
    Write a temp CSV file to be used as source catalog, similar to what
    is produced by the previous pipeline step, source detection.

    Parameters
    ----------
    tmp_path : pathlib.PosixPath
        A path-like object representing the path where to save the file.
    output_filename : string
        The output filename (with extension).
    catalog_data : numpy.ndarray, optional
        A numpy array with the (x, y) coordinates of the
        "detected" sources, by default None
    """
    header = ["x", "y"]
    # add shift
    src_detector_coords = catalog_data
    output = os.path.join(tmp_path, output_filename)
    with open(output, "w", encoding="UTF8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(src_detector_coords)


def add_tweakreg_catalog_attribute(tmp_path, input_dm, catalog_data=None):
    """
    Add tweakreg_catalog attribute to the meta, which is a mandatory
    attribute for TweakReg. The tweakreg_catalog attribute contains
    a table with the coordinates of the sources detected by the
    previous step. The coordinates will be used in the process of
    correcting the original WCS.

    Parameters
    ----------
    tmp_path : pathlib.PosixPath
        A path-like object representing the path where to save the file.
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.
    catalog_data : numpy.ndarray, optional
        A numpy array with the (x, y) coordinates of the
        "detected" sources, by default None (see note below).

    Note
    ----
    If no catalog_data is provided, a default catalog will be created
    by fetching data from Gaia within a search radius of 100 arcsec
    centered at RA=270, Dec=66.
    """
    tweakreg_catalog_filename = "base_image_sources.csv"
    if catalog_data is None:
        gaia_cat = get_catalog(ra=270, dec=66, sr=100 / 3600)
        gaia_source_coords = [
            (ra, dec) for ra, dec in zip(gaia_cat["ra"], gaia_cat["dec"])
        ]
        catalog_data = np.array(
            [
                input_dm.meta.wcs.world_to_pixel(ra, dec)
                for ra, dec in gaia_source_coords
            ]
        )
    create_base_image_source_catalog(
        tmp_path, tweakreg_catalog_filename, catalog_data=catalog_data
    )
    input_dm.meta["tweakreg_catalog"] = os.path.join(
        tmp_path, tweakreg_catalog_filename
    )


@pytest.fixture
def base_image():
    """
    Create a base image with a realistic WCSInfo and a WCS.

    Note
    ----
    The size of the image needs to be relatively large in order for
    the source catalog step to find a reasonable number of sources in the image.
    """

    def _base_image(shift_1=0, shift_2=0):
        l2 = maker_utils.mk_level2_image(shape=(2000, 2000))
        # update wcsinfo
        update_wcsinfo(l2)
        # add a dummy WCS object
        create_wcs_for_tweakreg_pipeline(l2, shift_1=shift_1, shift_2=shift_2)
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
    """
    Test that TweakReg raises an AttributeError if meta.tweakreg_catalog is missing.
    """
    img = base_image()
    with pytest.raises(Exception) as exec_info:
        TweakRegStep.call([img])

    assert type(exec_info.value) == AttributeError


def test_tweakreg_returns_modelcontainer(tmp_path, base_image):
    """Test that TweakReg returns a ModelContainer."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)
    res = TweakRegStep.call([img])

    assert type(res) == rdm.ModelContainer


def test_tweakreg_updates_cal_step(tmp_path, base_image):
    """Test that TweakReg updates meta.cal_step with tweakreg = COMPLETE."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)
    res = TweakRegStep.call([img])

    assert hasattr(res[0].meta.cal_step, "tweakreg")
    assert res[0].meta.cal_step.tweakreg == "COMPLETE"


def test_tweakreg_updates_group_id(tmp_path, base_image):
    """Test that TweakReg updates 'group_id' with a non-zero length string."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)
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
    """
    Test that TweakReg corrections are within tolerance.
    All the parametrized values are in arcsec.
    """
    img1 = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)
    img2 = request.getfixturevalue("base_image")(
        shift_1=1000 + shift_1 / 0.1, shift_2=1000 + shift_2 / 0.1
    )
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
    img = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)
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
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)

    step = TweakRegStep()
    step.save_abs_catalog = True
    step.abs_refcat = "my_ref_cat"

    with pytest.raises(Exception) as exec_info:
        step.process([img])

    assert type(exec_info.value) == ValueError


def test_tweakreg_compute_radius_and_fiducial():
    """
    Test that compute_radius returns correct values for footprint radius and fiducial.
    """
    # create a "header" of an image that's 100x100 pixels (@0.1"/pix in both axis)
    # and reference coordinates (ra=10 deg, dec=0 deg)
    # sitting at the origin of the pixel coordinates frame (0,0).
    img_shape = (100, 100)
    ref_val = (10, 0) * u.degree
    ref_pix = (0, 0)
    pix_scale = (0.1, 0.1) * u.arcsec
    theta = 0 * u.degree

    wcsobj = create_basic_wcs(
        img_shape=img_shape,
        ref_val=ref_val,
        ref_pix=ref_pix,
        pix_scale=pix_scale,
        theta=theta,
    )

    computed_radius, computed_fiducial = compute_radius(wcsobj)

    # expected radius is the euclidean distance
    # from the center to the furthest edge of the WCS
    ref_initial_coord = (
        np.array(
            [
                (wcsobj.bounding_box[0][0]),
                (wcsobj.bounding_box[1][0]),
            ]
        )
        * pix_scale.to("deg")[0].value
    )
    fiducial_coord = (
        np.array(
            [
                (wcsobj.bounding_box[0][0] + wcsobj.bounding_box[0][1]) / 2,
                (wcsobj.bounding_box[1][0] + wcsobj.bounding_box[1][1]) / 2,
            ]
        )
        * pix_scale.to("deg")[0].value
    )
    expected_radius = np.linalg.norm(fiducial_coord - ref_initial_coord)
    # expected fiducial is the center of the WCS
    expected_fiducial = wcsobj(img_shape[0] / 2, img_shape[1] / 2)

    assert_allclose(expected_radius, computed_radius)
    assert_allclose(expected_fiducial, computed_fiducial)


@pytest.mark.parametrize(
    "theta, offset_x, offset_y",
    [
        (0 * u.deg, 0, 0),
        (0 * u.deg, 1, 0),
        (0 * u.deg, 0, 1),
        (0 * u.deg, 1, 1),
        (0.05 * u.deg, 0, 0),
        (0.05 * u.deg, 1, 0),
        (0.05 * u.deg, 0, 1),
        (0.05 * u.deg, 1, 1),
        (0.1 * u.deg, 0, 0),
        (0.1 * u.deg, 1, 0),
        (0.1 * u.deg, 0, 1),
        (0.1 * u.deg, 1, 1),
    ],
)
def test_tweakreg_rotated_plane(tmp_path, theta, offset_x, offset_y, request):
    """
    Test that TweakReg returns accurate results.
    """
    gaia_cat = get_catalog(ra=270, dec=66, sr=100 / 3600)
    gaia_source_coords = [(ra, dec) for ra, dec in zip(gaia_cat["ra"], gaia_cat["dec"])]

    img = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)
    original_wcs = copy.deepcopy(img.meta.wcs)

    # calculate original (x,y) for Gaia sources
    original_xy_gaia_sources = np.array(
        [original_wcs.world_to_pixel(ra, dec) for ra, dec in gaia_source_coords]
    )
    # move Gaia sources around by applying linear transformations
    # to their coords in the projected plane (same as a "wrong WCS")
    rot_matrix = models.Rotation2D(angle=theta)
    transformed_xy_gaia_sources = np.array(
        [
            rot_matrix(x, y) + np.array([offset_x, offset_y])
            for x, y in original_xy_gaia_sources
        ]
    )
    # save modified catalog to meta.tweakreg_catalog
    # (coords in the projection plane)
    add_tweakreg_catalog_attribute(
        tmp_path, img, catalog_data=transformed_xy_gaia_sources
    )

    step = TweakRegStep()
    step.abs_minobj = 3
    step.process([img])

    # get world coords for Gaia sources using "wrong WCS"
    original_ref_source = [
        original_wcs.pixel_to_world(x, y) for x, y in transformed_xy_gaia_sources
    ]
    # get world coords for Gaia sources using tweaked WCS
    new_ref_source = [
        img.meta.wcs.pixel_to_world(x, y) for x, y in transformed_xy_gaia_sources
    ]
    # celestial coordinates for Gaia sources
    gaia_ref_source = [
        coord.SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
        for ra, dec in gaia_source_coords
    ]
    # calculate distance between "wrong WCS" result and Gaia
    # (rounded to the 10th decimal place to avoid floating point issues)
    dist1 = [
        np.round(gref.separation(oref), 10)
        for gref, oref in zip(gaia_ref_source, original_ref_source)
    ]
    # calculate distance between tweaked WCS result and Gaia
    # (rounded to the 10th decimal place to avoid floating point issues)
    dist2 = [
        np.round(gref.separation(nref), 10)
        for gref, nref in zip(gaia_ref_source, new_ref_source)
    ]

    assert np.array([np.less_equal(d2, d1) for d1, d2 in zip(dist1, dist2)]).all()
