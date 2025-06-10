import copy
import json
import os
import shutil
from io import StringIO

import numpy as np
import pytest
import requests
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from astropy.table import Table
from astropy.time import Time
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from numpy.random import default_rng
from roman_datamodels import datamodels as rdm
from roman_datamodels import stnode
from stcal.tweakreg.astrometric_utils import get_catalog

from romancal.datamodels import ModelLibrary
from romancal.tweakreg import tweakreg_step as trs
from romancal.tweakreg.tweakreg_step import _validate_catalog_columns


class MockConnectionError:
    def __init__(self, *args, **kwargs):
        raise requests.exceptions.ConnectionError


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
        t = Table(x.get("cat_data"), names=("x", "y"))
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
    ref_val: tuple[u.Quantity, u.Quantity] = (
        u.Quantity("10 deg"),
        u.Quantity("0 deg"),
    ),
    pix_scale: tuple[u.Quantity, u.Quantity] = (
        u.Quantity("0.1 arcsec"),
        u.Quantity("0.1 arcsec"),
    ),
    theta: u.Quantity | None = None,
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
    theta = u.Quantity("0 deg") if theta is None else theta

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
    wcs_obj.bounding_box = ((-0.5, shape[-2] + 0.5), (-0.5, shape[-1] + 0.5))

    input_dm.meta["wcs"] = wcs_obj


def get_catalog_data(input_dm, **kwargs):
    ra = kwargs.get("ra", 270)
    dec = kwargs.get("dec", 66)
    sr = kwargs.get("sr", 100 / 3600)
    add_shifts = kwargs.get("add_shifts", False)
    gaia_cat = get_catalog(right_ascension=ra, declination=dec, search_radius=sr)
    gaia_source_coords = [
        (ra, dec) for ra, dec in zip(gaia_cat["ra"], gaia_cat["dec"], strict=False)
    ]
    catalog_data = np.array(
        [
            input_dm.meta.wcs.world_to_pixel_values(ra, dec)
            for ra, dec in gaia_source_coords
        ]
    )
    if add_shifts:
        rng = np.random.default_rng(seed=int(ra + dec))
        shifts = rng.uniform(-1, 1, size=catalog_data.shape)
        catalog_data += shifts
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
    t = Table(src_detector_coords, names=("x", "y"))
    if save_catalogs:
        t.write((tmp_path / output), format=catalog_format)
    # mimic the same output format from SourceCatalogStep
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
    **kwargs,
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

    Notes
    ----
    - kwargs will be passed on to get_catalog_data();
    - if no catalog_data is provided, a default catalog will be created
    by fetching data from Gaia within a search radius of 100 arcsec
    centered at RA=270, Dec=66.
    """
    tweakreg_catalog_filename = catalog_filename
    if catalog_data is None:
        catalog_data = get_catalog_data(input_dm, **kwargs)

    source_catalog = create_base_image_source_catalog(
        tmp_path,
        tweakreg_catalog_filename,
        catalog_data=catalog_data,
        catalog_format=catalog_format,
        save_catalogs=save_catalogs,
    )

    input_dm.meta["source_catalog"] = stnode.SourceCatalog.create_fake_data()

    if save_catalogs:
        # SourceCatalogStep adds the catalog path+filename
        input_dm.meta.source_catalog["tweakreg_catalog_name"] = os.path.join(
            tmp_path, tweakreg_catalog_filename
        )
    else:
        # SourceCatalogStep attaches the catalog data as a structured array
        input_dm.meta.source_catalog["tweakreg_catalog"] = source_catalog


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
        l2 = rdm.ImageModel.create_fake_data(shape=(2000, 2000))
        l2.meta.filename = "none"
        l2.meta.cal_step = stnode.L2CalStep.create_fake_data()
        l2.meta.cal_logs = stnode.CalLogs.create_fake_data()
        l2.meta.exposure.start_time = Time("2016-01-01T00:00:00")
        # update wcsinfo
        update_wcsinfo(l2)
        # add a dummy WCS object
        create_wcs_for_tweakreg_pipeline(l2, shift_1=shift_1, shift_2=shift_2)
        return l2

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
    # make sure tweakreg_catalog_name doesn't exist
    img.meta.source_catalog = stnode.SourceCatalog.create_fake_data()
    assert "tweakreg_catalog_name" not in img.meta.source_catalog
    with pytest.raises(AttributeError):
        trs.TweakRegStep.call([img])


def test_tweakreg_returns_modellibrary_on_roman_datamodel_as_input(
    tmp_path, base_image
):
    """Test that TweakReg always returns a ModelLibrary when processing an open Roman DataModel as input."""

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename="img_1")

    test_input = img

    res = trs.TweakRegStep.call(test_input, save_l1_wcs=False)
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

    res = trs.TweakRegStep.call(test_input, save_l1_wcs=False)
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

    res = trs.TweakRegStep.call(test_input, save_l1_wcs=False)
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

    res = trs.TweakRegStep.call(test_input, save_l1_wcs=False)
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

    res = trs.TweakRegStep.call(test_input, save_l1_wcs=False)
    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert model.meta.cal_step.tweakreg == "COMPLETE"
            res.shelve(model, i, modify=False)


def test_tweakreg_updates_cal_step(tmp_path, base_image):
    """Test that TweakReg updates meta.cal_step with tweakreg = COMPLETE."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)
    res = trs.TweakRegStep.call([img], save_l1_wcs=False)

    with res:
        model = res.borrow(0)
        assert hasattr(model.meta.cal_step, "tweakreg")
        assert model.meta.cal_step.tweakreg == "COMPLETE"
        res.shelve(model, 0, modify=False)


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
        [img],
        save_abs_catalog=True,
        abs_refcat=abs_refcat,
        catalog_path=str(tmp_path),
        output_dir=str(tmp_path),
    )

    assert os.path.exists(tmp_path / abs_refcat_filename)


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
        [img],
        save_abs_catalog=True,
        abs_refcat=abs_refcat,
        catalog_path=str(tmp_path),
        output_dir=str(tmp_path),
    )

    assert os.path.exists(tmp_path / abs_refcat_filename)


def test_tweakreg_raises_error_on_invalid_abs_refcat(tmp_path, base_image):
    """Test that TweakReg raises an error when an invalid abs_refcat is provided."""
    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)

    with pytest.raises(TypeError):
        trs.TweakRegStep.call(
            [img],
            save_abs_catalog=True,
            abs_refcat="my_ref_cat",
            catalog_path=str(tmp_path),
            output_dir=str(tmp_path),
        )


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

    res = trs.TweakRegStep.call(
        asn_filepath,
        use_custom_catalogs=True,
        catalog_format=catalog_format,
        catfile=catfile,
        save_l1_wcs=False,
    )

    assert isinstance(res, ModelLibrary)

    with res:
        for i, (model, target) in enumerate(zip(res, [img1, img2, img3], strict=False)):
            assert hasattr(model.meta, "asn")

            assert model.meta.filename == target.meta.filename

            assert type(model) is type(target)

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
    gaia_cat = get_catalog(
        right_ascension=270, declination=66, search_radius=100 / 3600
    )
    gaia_source_coords = [
        (ra, dec) for ra, dec in zip(gaia_cat["ra"], gaia_cat["dec"], strict=False)
    ]

    img = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)
    original_wcs = copy.deepcopy(img.meta.wcs)

    # calculate original (x,y) for Gaia sources
    original_xy_gaia_sources = np.array(
        [original_wcs.world_to_pixel_values(ra, dec) for ra, dec in gaia_source_coords]
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

    trs.TweakRegStep.call([img], abs_minobj=3, save_l1_wcs=False)

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
        for gref, oref in zip(gaia_ref_source, original_ref_source, strict=False)
    ]
    # calculate distance between tweaked WCS result and Gaia
    # (rounded to the 10th decimal place to avoid floating point issues)
    dist2 = [
        np.round(gref.separation(nref), 10)
        for gref, nref in zip(gaia_ref_source, new_ref_source, strict=False)
    ]

    assert np.array(
        [np.less_equal(d2, d1) for d1, d2 in zip(dist1, dist2, strict=False)]
    ).all()


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

    res = trs.TweakRegStep.call(asn_filepath, save_l1_wcs=False)

    assert isinstance(res, ModelLibrary)

    with res:
        models = list(res)
        assert hasattr(models[0].meta, "asn")
        assert models[0].meta.asn["pool_name"] == asn_content["asn_pool"]
        assert models[1].meta.asn["pool_name"] == asn_content["asn_pool"]

        assert models[0].meta.filename == img_1.meta.filename
        assert models[1].meta.filename == img_2.meta.filename

        assert type(models[0]) is type(img_1)
        assert type(models[1]) is type(img_2)

        assert (models[0].data == img_1.data).all()
        assert (models[1].data == img_2.data).all()

        [res.shelve(m, i, modify=False) for i, m in enumerate(models)]


def test_fit_results_in_meta(tmp_path, base_image):
    """
    Test that the WCS fit results from tweakwcs are available in the meta tree.
    """

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img)

    res = trs.TweakRegStep.call([img], save_l1_wcs=False)

    assert isinstance(res, ModelLibrary)
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

    img1.meta.observation.program = 1
    img1.meta.observation["observation_id"] = "1"
    img2.meta.observation.program = 2
    img2.meta.observation["observation_id"] = "2"

    img1.meta["filename"] = "file1.asdf"
    img2.meta["filename"] = "file2.asdf"

    res = trs.TweakRegStep.call([img1, img2], save_l1_wcs=False)

    assert len(res.group_names) == 2


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
        x.meta.filename == y
        for x, y in zip(res_dict.get("datamodels"), catdict.keys(), strict=False)
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

    with pytest.raises(ValueError):
        trs._parse_catfile(catfile)


def test_update_source_catalog_coordinates(tmp_path, base_image):
    """Test that TweakReg updates the catalog coordinates with the tweaked WCS."""

    os.chdir(tmp_path)

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename="img_1")

    # create ImageSourceCatalogModel
    source_catalog = setup_source_catalog(img)
    source_catalog.write("img_1_cat.parquet", overwrite=True)

    # update tweakreg catalog name
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # run TweakRegStep
    res = trs.TweakRegStep.call([img])

    # tweak the current WCS using TweakRegStep and save the updated cat file
    with res:
        dm = res.borrow(0)
        assert dm.meta.source_catalog.tweakreg_catalog_name == "img_1_cat.parquet"
        update_catalog_coordinates(
            dm.meta.source_catalog.tweakreg_catalog_name, dm.meta.wcs
        )
        res.shelve(dm, 0)

    # read in saved catalog coords
    cat = Table.read("img_1_cat.parquet")
    cat_ra_centroid = cat["ra_centroid"]
    cat_dec_centroid = cat["dec_centroid"]
    cat_ra_psf = cat["ra_psf"]
    cat_dec_psf = cat["dec_psf"]

    # calculate world coords using tweaked WCS
    expected_centroid = img.meta.wcs(cat["xcentroid"], cat["ycentroid"])
    expected_psf = img.meta.wcs(cat["x_psf"], cat["y_psf"])

    # compare coordinates (make sure tweaked WCS was applied to cat file coords)
    np.testing.assert_array_equal(cat_ra_centroid, expected_centroid[0])
    np.testing.assert_array_equal(cat_dec_centroid, expected_centroid[1])
    np.testing.assert_array_equal(cat_ra_psf, expected_psf[0])
    np.testing.assert_array_equal(cat_dec_psf, expected_psf[1])


def test_source_catalog_coordinates_have_changed(tmp_path, base_image):
    """Test that the original catalog file content is different from the updated file."""

    os.chdir(tmp_path)

    img = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename="img_1")

    # create ImageSourceCatalogModel
    source_catalog = setup_source_catalog(img)
    source_catalog.write("img_1_cat.parquet", overwrite=True)

    # save original data
    shutil.copy("img_1_cat.parquet", "img_1_cat_original.parquet")

    # update tweakreg catalog name
    img.meta.source_catalog.tweakreg_catalog_name = "img_1_cat.parquet"

    # run TweakRegStep
    res = trs.TweakRegStep.call([img])

    # tweak the current WCS using TweakRegStep and save the updated cat file
    with res:
        dm = res.borrow(0)
        assert dm.meta.source_catalog.tweakreg_catalog_name == "img_1_cat.parquet"
        update_catalog_coordinates(
            dm.meta.source_catalog.tweakreg_catalog_name, dm.meta.wcs
        )
        res.shelve(dm, 0)

    cat_original = Table.read("img_1_cat_original.parquet")
    cat_updated = Table.read("img_1_cat.parquet")

    # set max absolute and relative tolerance to ~ 1/2 a pixel
    atol = u.Quantity(0.11 / 2, "arcsec").to("deg").value
    rtol = 5e-8

    # testing that nothing moved by more than 1/2 a pixel
    assert np.allclose(
        cat_original["ra_centroid"],
        cat_updated["ra_centroid"],
        atol=atol,
        rtol=rtol,
    )
    assert np.allclose(
        cat_original["dec_centroid"],
        cat_updated["dec_centroid"],
        atol=atol,
        rtol=rtol,
    )
    assert np.allclose(
        cat_original["ra_psf"],
        cat_updated["ra_psf"],
        atol=atol,
        rtol=rtol,
    )
    assert np.allclose(
        cat_original["dec_psf"],
        cat_updated["dec_psf"],
        atol=atol,
        rtol=rtol,
    )
    # testing that things did move by more than ~ 1/100 of a pixel
    assert not np.allclose(
        cat_original["ra_centroid"],
        cat_updated["ra_centroid"],
        atol=atol / 100,
        rtol=rtol / 100,
    )
    assert not np.allclose(
        cat_original["dec_centroid"],
        cat_updated["dec_centroid"],
        atol=atol / 100,
        rtol=rtol / 100,
    )
    assert not np.allclose(
        cat_original["ra_psf"],
        cat_updated["ra_psf"],
        atol=atol / 100,
        rtol=rtol / 100,
    )
    assert not np.allclose(
        cat_original["dec_psf"],
        cat_updated["dec_psf"],
        atol=atol / 100,
        rtol=rtol / 100,
    )


def setup_source_catalog(img):
    """
    Set up the source catalog.

    Notes
    -----
    This function reads the source catalog from a file, renames columns to match
    expected names, adds mock PSF coordinates, applies random shifts to the centroid
    and PSF coordinates, and calculates the world coordinates for the centroids.
    """
    # read in the mock table
    source_catalog = Table.read("img_1", format="ascii.ecsv")
    # rename columns to match expected column names
    source_catalog.rename_columns(["x", "y"], ["xcentroid", "ycentroid"])
    # add mock PSF coordinates
    source_catalog["x_psf"] = source_catalog["xcentroid"]
    source_catalog["y_psf"] = source_catalog["ycentroid"]

    # generate a set of random shifts to be added to the original coordinates
    seed = 13
    rng = default_rng(seed)
    shift_x = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    shift_y = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    # add random fraction of a pixel shifts to the centroid coordinates
    source_catalog["xcentroid"] += shift_x
    source_catalog["ycentroid"] += shift_y

    # generate another set of random shifts to be added to the original coordinates
    seed = 5
    rng = default_rng(seed)
    shift_x = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    shift_y = rng.uniform(-0.5, 0.5, size=len(source_catalog))
    # add random fraction of a pixel shifts to the centroid coordinates
    source_catalog["x_psf"] += shift_x
    source_catalog["y_psf"] += shift_y

    # calculate centroid world coordinates
    centroid = img.meta.wcs(
        source_catalog["xcentroid"],
        source_catalog["ycentroid"],
    )
    # calculate PSF world coordinates
    psf = img.meta.wcs(
        source_catalog["x_psf"],
        source_catalog["y_psf"],
    )
    # add world coordinates to catalog
    source_catalog["ra_centroid"], source_catalog["dec_centroid"] = centroid
    source_catalog["ra_psf"], source_catalog["dec_psf"] = psf
    # add units
    source_catalog["ra_centroid"].unit = u.deg
    source_catalog["dec_centroid"].unit = u.deg
    source_catalog["ra_psf"].unit = u.deg
    source_catalog["dec_psf"].unit = u.deg

    return source_catalog


def update_catalog_coordinates(tweakreg_catalog_name, tweaked_wcs):
    """
    Update the source catalog coordinates using the tweaked WCS.

    Parameters
    ----------
    tweakreg_catalog_name : str
        The name of the TweakReg catalog file produced by `SourceCatalog`.
    tweaked_wcs : `gwcs.wcs.WCS`
        The tweaked World Coordinate System (WCS) object.

    Returns
    -------
    None
    """
    # read in cat file
    catalog = Table.read(tweakreg_catalog_name)

    # define mapping between pixel and world coordinates
    colname_mapping = {
        ("xcentroid", "ycentroid"): ("ra_centroid", "dec_centroid"),
        ("x_psf", "y_psf"): ("ra_psf", "dec_psf"),
    }

    for k, v in colname_mapping.items():
        # get column names
        x_colname, y_colname = k
        ra_colname, dec_colname = v

        # calculate new coordinates using tweaked WCS and update catalog coordinates
        catalog[ra_colname], catalog[dec_colname] = tweaked_wcs(
            catalog[x_colname], catalog[y_colname]
        )

    catalog.write(tweakreg_catalog_name, overwrite=True)


@pytest.mark.parametrize(
    "exposure_type",
    ["WFI_GRISM", "WFI_PRISM", "WFI_DARK", "WFI_FLAT", "WFI_WFSC"],
)
def test_tweakreg_skips_invalid_exposure_types(exposure_type, tmp_path, base_image):
    """Test that TweakReg updates meta.cal_step with tweakreg = COMPLETE."""
    img1 = base_image(shift_1=1000, shift_2=1000)
    img1.meta.exposure.type = exposure_type
    img2 = base_image(shift_1=1000, shift_2=1000)
    img2.meta.exposure.type = exposure_type
    res = trs.TweakRegStep.call([img1, img2])

    assert isinstance(res, ModelLibrary)
    with res:
        for i, model in enumerate(res):
            assert hasattr(model.meta.cal_step, "tweakreg")
            assert model.meta.cal_step.tweakreg == "SKIPPED"
            res.shelve(model, i, modify=False)


@pytest.mark.parametrize(
    "catalog_data, expected_colnames, raises_exception",
    [
        # both 'x' and 'y' columns present
        ({"x": [1, 2, 3], "y": [4, 5, 6]}, ["x", "y"], False),
        # 'x_psf' and 'y_psf' columns present, should be renamed
        ({"x_psf": [1, 2, 3], "y_psf": [4, 5, 6]}, ["x", "y"], False),
        # 'x' present, 'y_psf' present, should rename 'y_psf' to 'y'
        ({"x": [1, 2, 3], "y_psf": [4, 5, 6]}, ["x", "y"], False),
        # 'x_psf' present, 'y' present, should rename 'x_psf' to 'x'
        ({"x_psf": [1, 2, 3], "y": [4, 5, 6]}, ["x", "y"], False),
        # neither 'x' nor 'x_psf' present
        ({"y": [4, 5, 6]}, None, True),
        # neither 'y' nor 'y_psf' present
        ({"x": [1, 2, 3]}, None, True),
        # no relevant columns present
        (
            {"a": [1, 2, 3], "b": [4, 5, 6]},
            None,
            True,
        ),
    ],
)
def test_validate_catalog_columns(catalog_data, expected_colnames, raises_exception):
    """Test that TweakRegStep._validate_catalog_columns() correctly validates the
    presence of required columns ('x' and 'y') in the provided catalog."""
    catalog = Table(catalog_data)
    if raises_exception:
        with pytest.raises(ValueError):
            _validate_catalog_columns(catalog)
    else:
        _validate_catalog_columns(catalog)
        assert set(catalog.colnames) == set(expected_colnames)


def test_tweakreg_handles_mixed_exposure_types(tmp_path, base_image):
    """Test that TweakReg can handle mixed exposure types
    (non-WFI_IMAGE data will be marked as SKIPPED only and won't be processed)."""
    img1 = base_image(shift_1=1000, shift_2=1000)
    img1.meta.exposure.type = "WFI_GRISM"

    img2 = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img2, catalog_filename="img2")
    img2.meta.exposure.type = "WFI_IMAGE"

    img3 = base_image(shift_1=1000, shift_2=1000)
    img3.meta.exposure.type = "WFI_GRISM"

    img4 = base_image(shift_1=1000, shift_2=1000)
    add_tweakreg_catalog_attribute(tmp_path, img4, catalog_filename="img4")
    img4.meta.exposure.type = "WFI_IMAGE"

    img5 = base_image(shift_1=1000, shift_2=1000)
    img5.meta.exposure.type = "WFI_GRISM"

    res = trs.TweakRegStep.call([img1, img2, img3, img4, img5])

    assert len(res) == 5
    assert img1.meta.cal_step.tweakreg == "SKIPPED"
    assert img2.meta.cal_step.tweakreg == "COMPLETE"
    assert img3.meta.cal_step.tweakreg == "SKIPPED"
    assert img4.meta.cal_step.tweakreg == "COMPLETE"
    assert img5.meta.cal_step.tweakreg == "SKIPPED"


def test_tweakreg_updates_s_region(tmp_path, base_image):
    """Test that the TweakRegStep updates the s_region attribute."""
    img = base_image(shift_1=1000, shift_2=1000)
    old_fake_s_region = "POLYGON ICRS 1.0000000000000 2.0000000000000 3.0000000000000 4.0000000000000 5.0000000000000 6.0000000000000 7.0000000000000 8.0000000000000 "
    img.meta.wcsinfo["s_region"] = old_fake_s_region
    add_tweakreg_catalog_attribute(tmp_path, img, catalog_filename="img")

    # call TweakRegStep to update WCS & S_REGION
    res = trs.TweakRegStep.call([img])

    with res:
        for i, model in enumerate(res):
            assert model.meta.wcsinfo.s_region != old_fake_s_region
            res.shelve(model, i, modify=False)
