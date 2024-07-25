import copy
import json
import os
from io import StringIO
from pathlib import Path
from typing import Tuple

import numpy as np
import pytest
import requests
from astropy import coordinates as coord
from astropy import table
from astropy import units as u
from astropy.modeling import models
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from astropy.table import Table
from astropy.time import Time
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.datamodels import ModelLibrary
from romancal.tweakreg import tweakreg_step as trs
from romancal.tweakreg.astrometric_utils import get_catalog


class MockConnectionError:
    def __init__(self, *args, **kwargs):
        raise requests.exceptions.ConnectionError


@pytest.fixture()
def test_data_dir():
    return Path.joinpath(Path(__file__).parent, "data")


def create_custom_catalogs(tmp_path, base_image, catalog_format="ascii.ecsv"):
    """
    Creates a custom catalog with three datamodels and its respective catalog data.
    """
    img1 = base_image(shift_1=1000, shift_2=1000)
    img2 = base_image(shift_1=1010, shift_2=1010)
    img3 = base_image(shift_1=1020, shift_2=1020)

    img1.meta.filename = "img1.asdf"
    img2.meta.filename = "img2.asdf"
    img3.meta.filename = "img3.asdf"

    # create valid custom catalog data to be used with each input datamodel
    catalog_data1 = get_catalog_data(img1)
    catalog_data2 = get_catalog_data(img2)
    catalog_data3 = get_catalog_data(img3)

    custom_catalog_map = [
        {
            "cat_filename": "ref_catalog_1",
            "cat_datamodel": img1.meta.filename,
            "cat_data": catalog_data1,
        },
        {
            "cat_filename": "ref_catalog_2",
            "cat_datamodel": img2.meta.filename,
            "cat_data": catalog_data2,
        },
        {
            "cat_filename": "ref_catalog_3",
            "cat_datamodel": img3.meta.filename,
            "cat_data": catalog_data3,
        },
    ]

    # create catfile
    catfile = str(tmp_path / "catfile.txt")
    catfile_content = StringIO()
    for x in custom_catalog_map:
        # write line to catfile
        catfile_content.write(f"{x.get('cat_datamodel')} {x.get('cat_filename')}\n")
        # write out the catalog data
        t = table.Table(x.get("cat_data"), names=("x", "y"))
        t.write(tmp_path / x.get("cat_filename"), format=catalog_format)
    with open(catfile, mode="w") as f:
        print(catfile_content.getvalue(), file=f)

    return {"catfile": catfile, "datamodels": [img1, img2, img3]}


def create_asn_file(tmp_path, members_mapping=None):
    asn_content = """
        {
            "asn_type": "None",
            "asn_rule": "DMS_ELPP_Base",
            "version_id": null,
            "code_version": "0.9.1.dev28+ge987cc9.d20230106",
            "degraded_status": "No known degraded exposures in association.",
            "program": "noprogram",
            "constraints": "No constraints",
            "asn_id": "a3001",
            "target": "none",
            "asn_pool": "test_pool_name",
            "products": [
                {
                    "name": "files.asdf",
                    "members": [
                        {
                            "expname": "img_1.asdf",
                            "exptype": "science"
                        },
                        {
                            "expname": "img_2.asdf",
                            "exptype": "science"
                        }
                    ]
                }
            ]
        }
    """
    if members_mapping is not None:
        asn_dict = json.loads(asn_content)
        asn_dict["products"][0]["members"] = []
        for x in members_mapping:
            asn_dict["products"][0]["members"].append(x)
        asn_content = json.dumps(asn_dict)

    asn_file_path = str(tmp_path / "sample_asn.json")
    asn_file = StringIO()
    asn_file.write(asn_content)
    with open(asn_file_path, mode="w") as f:
        print(asn_file.getvalue(), file=f)

    return asn_file_path


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
        "POLYGON ICRS "
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
    ref_val: Tuple[u.Quantity, u.Quantity] = (
        u.Quantity("10 deg"),
        u.Quantity("0 deg"),
    ),
    pix_scale: Tuple[u.Quantity, u.Quantity] = (
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
    wcsobj.bounding_box = (
        (-0.5, img_shape[0] + 0.5),
        (-0.5, img_shape[1] + 0.5),
    )

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


def get_catalog_data(input_dm):
    gaia_cat = get_catalog(ra=270, dec=66, sr=100 / 3600)
    gaia_source_coords = [(ra, dec) for ra, dec in zip(gaia_cat["ra"], gaia_cat["dec"])]
    catalog_data = np.array(
        [input_dm.meta.wcs.world_to_pixel(ra, dec) for ra, dec in gaia_source_coords]
    )
    return catalog_data


def create_base_image_source_catalog(
    tmp_path,
    output_filename,
    catalog_data,
    catalog_format: str = "ascii.ecsv",
    save_catalogs=True,
):
    """
    Write a temp CSV file to be used as source catalog, similar to what
    is produced by the previous pipeline step, source detection.

    Parameters
    ----------
    tmp_path : pathlib.PosixPath
        A path-like object representing the path where to save the file.
    output_filename : string
        The output filename (with extension).
    catalog_data : numpy.ndarray
        A numpy array with the (x, y) coordinates of the
        "detected" sources
    catalog_format : str, optional
        A string indicating the catalog format.
    save_catalogs : boolean, optional
        A boolean indicating whether the source catalog should be saved to disk.
    """
    src_detector_coords = catalog_data
    output = os.path.join(tmp_path, output_filename)
    t = table.Table(src_detector_coords, names=("x", "y"))
    if save_catalogs:
        t.write((tmp_path / output), format=catalog_format)
    # mimic the same output format from SourceDetectionStep
    t.add_column([i for i in range(len(t))], name="id", index=0)
    t.add_column([np.float64(i) for i in range(len(t))], name="flux")
    t.rename_columns(["x", "y"], ["xcentroid", "ycentroid"])
    return t.as_array()


def add_tweakreg_catalog_attribute(
    tmp_path,
    input_dm,
    catalog_filename="base_image_sources",
    catalog_data=None,
    catalog_format: str = "ascii.ecsv",
    save_catalogs=True,
):
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
    catalog_filename : str
        Filename to be used when saving the catalog to disk.
    catalog_data : numpy.ndarray, optional
        A numpy array with the (x, y) coordinates of the
        "detected" sources, by default None (see note below).
    catalog_format : str, optional
        A string indicating the catalog format.
    save_catalogs : boolean, optional
        A boolean indicating whether the source catalog should be saved to disk.

    Note
    ----
    If no catalog_data is provided, a default catalog will be created
    by fetching data from Gaia within a search radius of 100 arcsec
    centered at RA=270, Dec=66.
    """
    tweakreg_catalog_filename = catalog_filename
    if catalog_data is None:
        catalog_data = get_catalog_data(input_dm)

    source_catalog = create_base_image_source_catalog(
        tmp_path,
        tweakreg_catalog_filename,
        catalog_data=catalog_data,
        catalog_format=catalog_format,
        save_catalogs=save_catalogs,
    )

    input_dm.meta["source_detection"] = maker_utils.mk_source_detection()

    if save_catalogs:
        # SourceDetectionStep adds the catalog path+filename
        input_dm.meta.source_detection["tweakreg_catalog_name"] = os.path.join(
            tmp_path, tweakreg_catalog_filename
        )
    else:
        # SourceDetectionStep attaches the catalog data as a structured array
        input_dm.meta.source_detection["tweakreg_catalog"] = source_catalog


@pytest.fixture
def base_image():
    """
    Create a base image with a realistic WCSInfo and a WCS.

    Notes
    -----
    The size of the image needs to be relatively large in order for
    the source catalog step to find a reasonable number of sources in the image.

    shift_1 and shift_2 (units in pixel) are used to shift the WCS projection plane.
    """

    def _base_image(shift_1=0, shift_2=0):
        l2 = maker_utils.mk_level2_image(shape=(2000, 2000))
        l2.meta.exposure.mid_time = Time("2016-01-01T00:00:00")
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
        (list(), (Exception,)),
        ([""], (Exception,)),
        (["", ""], (Exception,)),
        ("", (Exception,)),
        ([1, 2, 3], (Exception,)),
    ],
)
def test_tweakreg_raises_error_on_invalid_input(input, error_type):
    # sourcery skip: list-literal
    """Test that TweakReg raises an error when an invalid input is provided."""
    with pytest.raises(error_type):
        trs.TweakRegStep.call(input)


def test_tweakreg_raises_attributeerror_on_missing_tweakreg_catalog(base_image):
    """
    Test that TweakReg raises an AttributeError if meta.tweakreg_catalog is missing.
    """
    img = base_image()
    with pytest.raises(AttributeError):
        trs.TweakRegStep.call([img])


def test_tweakreg_returns_modellibrary_on_roman_datamodel_as_input(
    tmp_path, base_image
):
    """Test that TweakReg always returns a ModelLibrary when processing an open Roman DataModel as input."""

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename="img_1")

    test_input = img

    res = trs.TweakRegStep.call(test_input)
    assert isinstance(res, ModelLibrary)
    with res:
        model = res.borrow(0)
        assert model.meta.cal_step.tweakreg == "COMPLETE"
        res.shelve(model, 0, modify=False)


def test_tweakreg_returns_modellibrary_on_modellibrary_as_input(tmp_path, base_image):
    """Test that TweakReg always returns a ModelLibrary when processing a ModelLibrary as input."""

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename="img_1")

    test_input = ModelLibrary([img])

    res = trs.TweakRegStep.call(test_input)
    assert isinstance(res, ModelLibrary)
    with res:
        model = res.borrow(0)
        assert model.meta.cal_step.tweakreg == "COMPLETE"
        res.shelve(model, 0, modify=False)


def test_tweakreg_returns_modellibrary_on_association_file_as_input(
    tmp_path, base_image
):
    """Test that TweakReg always returns a ModelLibrary when processing an association file as input."""

    img_1 = base_image(shift_1=1000, shift_2=1000)
    img_2 = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    add_tweakreg_catalog_attribute(tmp_path, img_2, catalog_filename="img_2")
    img_1.save(tmp_path / "img_1.asdf")
    img_2.save(tmp_path / "img_2.asdf")
    asn_filepath = create_asn_file(tmp_path)

    test_input = asn_filepath

    res = trs.TweakRegStep.call(test_input)
    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert model.meta.cal_step.tweakreg == "COMPLETE"
            res.shelve(model, i, modify=False)


def test_tweakreg_returns_modellibrary_on_list_of_asdf_file_as_input(
    tmp_path, base_image
):
    """Test that TweakReg always returns a ModelLibrary when processing a list of ASDF files as input."""

    img_1 = base_image(shift_1=1000, shift_2=1000)
    img_2 = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    add_tweakreg_catalog_attribute(tmp_path, img_2, catalog_filename="img_2")
    img_1.save(tmp_path / "img_1.asdf")
    img_2.save(tmp_path / "img_2.asdf")

    tmp_path_str = tmp_path.as_posix()
    test_input = [
        f"{tmp_path_str}/img_1.asdf",
        f"{tmp_path_str}/img_2.asdf",
    ]

    res = trs.TweakRegStep.call(test_input)
    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert model.meta.cal_step.tweakreg == "COMPLETE"
            res.shelve(model, i, modify=False)


def test_tweakreg_returns_modellibrary_on_list_of_roman_datamodels_as_input(
    tmp_path, base_image
):
    """Test that TweakReg always returns a ModelLibrary when processing a list of open Roman datamodels as input."""
    img_1 = base_image(shift_1=1000, shift_2=1000)
    img_2 = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    add_tweakreg_catalog_attribute(tmp_path, img_2, catalog_filename="img_2")
    img_1.save(tmp_path / "img_1.asdf")
    img_2.save(tmp_path / "img_2.asdf")

    test_input = [img_1, img_2]

    res = trs.TweakRegStep.call(test_input)
    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert model.meta.cal_step.tweakreg == "COMPLETE"
            res.shelve(model, i, modify=False)


def test_tweakreg_updates_cal_step(tmp_path, base_image):
    """Test that TweakReg updates meta.cal_step with tweakreg = COMPLETE."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)
    res = trs.TweakRegStep.call([img])

    with res:
        model = res.borrow(0)
        assert hasattr(model.meta.cal_step, "tweakreg")
        assert model.meta.cal_step.tweakreg == "COMPLETE"
        res.shelve(model, 0, modify=False)


def test_tweakreg_updates_group_id(tmp_path, base_image):
    """Test that TweakReg updates 'group_id' with a non-zero length string."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)
    res = trs.TweakRegStep.call([img])

    with res:
        model = res.borrow(0)
        assert hasattr(model.meta, "group_id")
        res.shelve(model, 0, modify=False)


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

    step = trs.TweakRegStep()
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

    res = trs._common_name(img_list)

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
        trs._common_name(filename_list)

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
    catalog_filename = "ref_catalog.ecsv"
    abs_refcat_filename = f"fit_{abs_refcat.lower()}_ref.ecsv"
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename=catalog_filename)

    trs.TweakRegStep.call(
        [img], save_abs_catalog=True, abs_refcat=abs_refcat, catalog_path=str(tmp_path)
    )

    assert os.path.exists(tmp_path / abs_refcat_filename)
    # clean up
    os.remove(tmp_path / abs_refcat_filename)
    os.remove(tmp_path / catalog_filename)


@pytest.mark.parametrize(
    "abs_refcat",
    (None, ""),
)
def test_tweakreg_defaults_to_valid_abs_refcat(tmp_path, abs_refcat, request):
    """Test that TweakReg defaults to DEFAULT_ABS_REFCAT on invalid values."""
    img = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)
    catalog_filename = "ref_catalog.ecsv"
    abs_refcat_filename = f"fit_{trs.DEFAULT_ABS_REFCAT.lower()}_ref.ecsv"
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename=catalog_filename)

    trs.TweakRegStep.call(
        [img], save_abs_catalog=True, abs_refcat=abs_refcat, catalog_path=str(tmp_path)
    )

    assert os.path.exists(tmp_path / abs_refcat_filename)
    # clean up
    os.remove(tmp_path / abs_refcat_filename)
    os.remove(tmp_path / catalog_filename)


def test_tweakreg_raises_error_on_invalid_abs_refcat(tmp_path, base_image):
    """Test that TweakReg raises an error when an invalid abs_refcat is provided."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)

    with pytest.raises(Exception) as exec_info:
        trs.TweakRegStep.call([img], save_abs_catalog=True, abs_refcat="my_ref_cat")

    assert type(exec_info.value) == ValueError


def test_tweakreg_combine_custom_catalogs_and_asn_file(tmp_path, base_image):
    """
    Test that TweakRegStep can handle a custom catalog for the members of an ASN file.
    In this case, the user can create a custom catalog file (catfile) for each of the
    members of an ASN file.
    """
    # create custom catalog and input datamodels
    catalog_format = "ascii.ecsv"
    res_dict = create_custom_catalogs(
        tmp_path, base_image, catalog_format=catalog_format
    )
    catfile = res_dict.get("catfile")
    img1, img2, img3 = res_dict.get("datamodels")
    add_tweakreg_catalog_attribute(tmp_path, img1, catalog_filename="img1")
    add_tweakreg_catalog_attribute(tmp_path, img2, catalog_filename="img2")
    add_tweakreg_catalog_attribute(tmp_path, img3, catalog_filename="img3")
    img1.save(tmp_path / "img1.asdf")
    img2.save(tmp_path / "img2.asdf")
    img3.save(tmp_path / "img3.asdf")

    # create ASN file
    asn_filepath = create_asn_file(
        tmp_path,
        members_mapping=[
            {"expname": img1.meta.filename, "exptype": "science"},
            {"expname": img2.meta.filename, "exptype": "science"},
            {"expname": img3.meta.filename, "exptype": "science"},
        ],
    )
    with open(asn_filepath) as f:
        asn_content = json.load(f)

    res = trs.TweakRegStep.call(
        asn_filepath,
        use_custom_catalogs=True,
        catalog_format=catalog_format,
        catfile=catfile,
    )

    assert type(res) == ModelLibrary

    with res:
        for i, (model, target) in enumerate(zip(res, [img1, img2, img3])):
            assert hasattr(model.meta, "asn")

            assert (
                model.meta["exptype"]
                == asn_content["products"][0]["members"][i]["exptype"]
            )

            assert model.meta.filename == target.meta.filename

            assert type(model) == type(target)

            assert (model.data == target.data).all()

            res.shelve(model, i, modify=False)


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

    trs.TweakRegStep.call([img], abs_minobj=3)

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


def test_tweakreg_parses_asn_correctly(tmp_path, base_image):
    """Test that TweakReg can parse an ASN file properly."""

    img_1 = base_image(shift_1=1000, shift_2=1000)
    img_2 = base_image(shift_1=1000, shift_2=1000)
    img_1.meta["filename"] = "img_1.asdf"
    img_2.meta["filename"] = "img_2.asdf"
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    add_tweakreg_catalog_attribute(tmp_path, img_2, catalog_filename="img_2")
    img_1.save(tmp_path / "img_1.asdf")
    img_2.save(tmp_path / "img_2.asdf")
    asn_filepath = create_asn_file(tmp_path)
    with open(asn_filepath) as f:
        asn_content = json.load(f)

    res = trs.TweakRegStep.call(asn_filepath)

    assert type(res) == ModelLibrary

    with res:
        models = list(res)
        assert hasattr(models[0].meta, "asn")
        assert (
            models[0].meta["exptype"]
            == asn_content["products"][0]["members"][0]["exptype"]
        )
        assert (
            models[1].meta["exptype"]
            == asn_content["products"][0]["members"][1]["exptype"]
        )
        assert models[0].meta.asn["pool_name"] == asn_content["asn_pool"]
        assert models[1].meta.asn["pool_name"] == asn_content["asn_pool"]

        assert models[0].meta.filename == img_1.meta.filename
        assert models[1].meta.filename == img_2.meta.filename

        assert type(models[0]) == type(img_1)
        assert type(models[1]) == type(img_2)

        assert (models[0].data == img_1.data).all()
        assert (models[1].data == img_2.data).all()

        [res.shelve(m, i, modify=False) for i, m in enumerate(models)]


def test_tweakreg_raises_error_on_connection_error_to_the_vo_service(
    tmp_path, base_image, monkeypatch
):
    """
    Test that TweakReg raises an error when there is a connection error with
    the VO API server, which means that an absolute reference catalog cannot be created.
    """

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)

    monkeypatch.setattr("requests.get", MockConnectionError)
    res = trs.TweakRegStep.call([img])

    assert type(res) == ModelLibrary
    assert len(res) == 1
    with res:
        model = res.borrow(0)
        assert model.meta.cal_step.tweakreg.lower() == "skipped"
        res.shelve(model, 0, modify=False)


def test_fit_results_in_meta(tmp_path, base_image):
    """
    Test that the WCS fit results from tweakwcs are available in the meta tree.
    """

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)

    res = trs.TweakRegStep.call([img])

    assert type(res) == ModelLibrary
    with res:
        for i, model in enumerate(res):
            assert hasattr(model.meta, "wcs_fit_results")
            assert len(model.meta.wcs_fit_results) > 0
            res.shelve(model, i, modify=False)


def test_tweakreg_handles_multiple_groups(tmp_path, base_image):
    """
    Test that TweakRegStep can perform relative alignment for all images in the groups
    before performing absolute alignment.
    """
    img1 = base_image(shift_1=1000, shift_2=1000)
    img2 = base_image(shift_1=990, shift_2=990)
    add_tweakreg_catalog_attribute(tmp_path, img1, catalog_filename="img1")
    add_tweakreg_catalog_attribute(tmp_path, img2, catalog_filename="img2")

    img1.meta.observation["program"] = "-program_id1"
    img2.meta.observation["program"] = "-program_id2"

    img1.meta["filename"] = "file1.asdf"
    img2.meta["filename"] = "file2.asdf"

    res = trs.TweakRegStep.call([img1, img2])

    assert len(res.group_names) == 2
    with res:
        for r, i in zip(res, [img1, img2]):
            assert (
                r.meta.group_id.split("-")[1]
                == i.meta.observation.program.split("-")[1]
            )
            res.shelve(r, modify=False)


@pytest.mark.parametrize(
    "column_names",
    [("x", "y"), ("xcentroid", "ycentroid")],
)
def test_imodel2wcsim_valid_column_names(tmp_path, base_image, column_names):
    """
    Test that _imodel2wcsim handles different catalog column names.
    """
    img_1 = base_image(shift_1=1000, shift_2=1000)
    img_2 = base_image(shift_1=1030, shift_2=1030)
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    add_tweakreg_catalog_attribute(tmp_path, img_2, catalog_filename="img_2")
    # set meta.tweakreg_catalog (this is automatically added by TweakRegStep)
    catalog_format = "ascii.ecsv"
    for x in [img_1, img_2]:
        x.meta["tweakreg_catalog"] = Table.read(
            x.meta.source_detection.tweakreg_catalog_name,
            format=catalog_format,
        )
        x.meta.tweakreg_catalog.rename_columns(("x", "y"), column_names)
    xname, yname = column_names

    images = ModelLibrary([img_1, img_2])

    step = trs.TweakRegStep()
    with images:
        for i, (m, target) in enumerate(zip(images, [img_1, img_2])):
            imcat = step._imodel2wcsim(m)
            assert (
                imcat.meta["catalog"]["x"] == target.meta.tweakreg_catalog[xname]
            ).all()
            assert (
                imcat.meta["catalog"]["y"] == target.meta.tweakreg_catalog[yname]
            ).all()
            images.shelve(m, i, modify=False)


@pytest.mark.parametrize(
    "column_names",
    [
        ("x_centroid", "y_centroid"),
        ("x_cen", "y_cen"),
    ],
)
def test_imodel2wcsim_error_invalid_column_names(tmp_path, base_image, column_names):
    """
    Test that _imodel2wcsim raises a ValueError on invalid catalog column names.
    """
    img_1 = base_image(shift_1=1000, shift_2=1000)
    img_2 = base_image(shift_1=1030, shift_2=1030)
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    add_tweakreg_catalog_attribute(tmp_path, img_2, catalog_filename="img_2")
    # set meta.tweakreg_catalog (this is automatically added by TweakRegStep)
    catalog_format = "ascii.ecsv"
    for x in [img_1, img_2]:
        x.meta["tweakreg_catalog"] = Table.read(
            x.meta.source_detection.tweakreg_catalog_name,
            format=catalog_format,
        )
        x.meta.tweakreg_catalog.rename_columns(("x", "y"), column_names)

    images = ModelLibrary([img_1, img_2])

    step = trs.TweakRegStep()
    with pytest.raises(ValueError):
        with images:
            for i, model in enumerate(images):
                images.shelve(model, i, modify=False)
                step._imodel2wcsim(model)


def test_imodel2wcsim_error_invalid_catalog(tmp_path, base_image):
    """
    Test that _imodel2wcsim raises an error on invalid catalog format.
    """
    img_1 = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img_1, catalog_filename="img_1")
    # set meta.tweakreg_catalog (this is automatically added by TweakRegStep)
    img_1.meta["tweakreg_catalog"] = "nonsense"

    images = ModelLibrary([img_1])

    step = trs.TweakRegStep()
    with pytest.raises(AttributeError):
        with images:
            for i, model in enumerate(images):
                images.shelve(model, i, modify=False)
                step._imodel2wcsim(model)


def test_parse_catfile_valid_catalog(tmp_path, base_image):
    """
    Test that _parse_catfile can parse a custom catalog with valid format.
    """
    # create custom catalog file and input datamodels
    catalog_format = "ascii.ecsv"
    res_dict = create_custom_catalogs(
        tmp_path, base_image, catalog_format=catalog_format
    )
    catfile = res_dict.get("catfile")
    catdict = trs._parse_catfile(catfile)

    assert all(
        x.meta.filename == y for x, y in zip(res_dict.get("datamodels"), catdict.keys())
    )


@pytest.mark.parametrize("catfile", (None, ""))
def test_parse_catfile_returns_none(catfile):
    """
    Test that _parse_catfile returns None when catfile = None or catfile = "".
    """
    catdict = trs._parse_catfile(catfile=catfile)

    assert catdict is None


@pytest.mark.parametrize(
    "catfile_line_content",
    ["img1.asdf\nimg2.asdf\nimg3.asdf"],
)
def test_parse_catfile_returns_none_on_invalid_content(tmp_path, catfile_line_content):
    """
    Test that _parse_catfile returns a dict where all the values are None
    if only filename is present in catfile (i.e. no associated catalog).
    """
    # create custom catalog file and input datamodels
    catfile = str(tmp_path / "catfile.txt")
    catfile_content = StringIO()
    # write empty line to catfile
    catfile_content.write(catfile_line_content)
    # write StringIO object to disk
    with open(catfile, mode="w") as f:
        print(catfile_content.getvalue(), file=f)

    catdict = trs._parse_catfile(catfile)

    assert not all(catdict.values())


@pytest.mark.parametrize(
    "catfile_line_content",
    ["img1.asdf column1 column2 column3"],
)
def test_parse_catfile_raises_error_on_invalid_content(tmp_path, catfile_line_content):
    """
    Test that _parse_catfile raises an error if catfile contains more than
    two columns (i.e. image name and corresponding catalog path).
    """
    # create custom catalog file and input datamodels
    catfile = str(tmp_path / "catfile.txt")
    catfile_content = StringIO()
    # write empty line to catfile
    catfile_content.write(catfile_line_content)
    # write StringIO object to disk
    with open(catfile, mode="w") as f:
        print(catfile_content.getvalue(), file=f)

    with pytest.raises(Exception) as exec_info:
        trs._parse_catfile(catfile)

    assert type(exec_info.value) == ValueError
