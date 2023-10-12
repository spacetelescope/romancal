import os
from io import StringIO
from typing import Tuple

import numpy as np
import pytest
import requests
from astropy import coordinates as coord
from astropy import table
from astropy import units as u
from astropy.modeling import models
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from astropy.stats import mad_std
from astropy.time import Time
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from numpy.testing import assert_allclose
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.tweakreg.astrometric_utils import (
    compute_radius,
    create_astrometric_catalog,
    get_catalog,
)

ARAD = np.pi / 180.0


class MockConnectionError:
    def __init__(self, *args, **kwargs):
        raise requests.exceptions.ConnectionError


def get_parallax_correction_barycenter(epoch, gaia_ref_epoch_coords):
    """
    Calculates the parallax correction in the Earth barycenter frame for a given epoch
    and Gaia reference epoch coordinates (i.e. Gaia coordinates at the reference epoch).

    Parameters
    ----------
    epoch : float
        The epoch for which the parallax correction is calculated.
    gaia_ref_epoch_coords : dict
        The Gaia reference epoch coordinates, including 'ra', 'dec', and 'parallax'.

    Returns
    -------
    tuple
        A tuple containing the delta_ra and delta_dec values of the parallax correction
        in degrees.

    Examples
    --------
    .. code-block :: python
        epoch = 2022.5
        gaia_coords = {'ra': 180.0, 'dec': 45.0, 'parallax': 10.0}
        correction = get_parallax_correction_earth_barycenter(epoch, gaia_coords)
        print(correction)
        (0.001, -0.002)
    """  # noqa: E501

    obs_date = Time(epoch, format="decimalyear")
    earths_center_barycentric_coords = coord.get_body_barycentric(
        "earth", obs_date, ephemeris="builtin"
    )
    earth_X = earths_center_barycentric_coords.x
    earth_Y = earths_center_barycentric_coords.y
    earth_Z = earths_center_barycentric_coords.z

    # angular displacement components
    # (see eq. 8.15 of "Spherical Astronomy" by Robert M. Green)
    delta_ra = (
        u.Quantity(gaia_ref_epoch_coords["parallax"], unit="mas").to(u.rad)
        * (1 / np.cos(gaia_ref_epoch_coords["dec"] * ARAD))
        * (
            earth_X.value * np.sin(gaia_ref_epoch_coords["ra"] * ARAD)
            - earth_Y.value * np.cos(gaia_ref_epoch_coords["ra"] * ARAD)
        )
    ).to("deg")
    delta_dec = (
        u.Quantity(gaia_ref_epoch_coords["parallax"], unit="mas").to(u.rad)
        * (
            earth_X.value
            * np.cos(gaia_ref_epoch_coords["ra"] * ARAD)
            * np.sin(gaia_ref_epoch_coords["dec"] * ARAD)
            + earth_Y.value
            * np.sin(gaia_ref_epoch_coords["ra"] * ARAD)
            * np.sin(gaia_ref_epoch_coords["dec"] * ARAD)
            - earth_Z.value * np.cos(gaia_ref_epoch_coords["dec"] * ARAD)
        )
    ).to("deg")

    return delta_ra, delta_dec


def get_proper_motion_correction(epoch, gaia_ref_epoch_coords, gaia_ref_epoch):
    """
    Calculates the proper motion correction for a given epoch and Gaia reference epoch
    coordinates.

    Parameters
    ----------
    epoch : float
        The epoch for which the proper motion correction is calculated.
    gaia_ref_epoch_coords : dict
        A dictionary containing Gaia reference epoch coordinates.
    gaia_ref_epoch : float
        The Gaia reference epoch.

    Returns
    -------
    None

    Examples
    --------
    .. code-block:: python
    epoch = 2022.5
    gaia_coords = {
        "ra": 180.0,
        "dec": 45.0,
        "pmra": 2.0,
        "pmdec": 1.5
    }
    gaia_ref_epoch = 2020.0
    get_proper_motion_correction(epoch, gaia_coords, gaia_ref_epoch)
    """  # noqa: E501

    expected_new_dec = (
        np.array(
            gaia_ref_epoch_coords["dec"] * 3600
            + (epoch - gaia_ref_epoch) * gaia_ref_epoch_coords["pmdec"] / 1000
        )
        / 3600
    )
    average_dec = np.array(
        [
            np.mean([new, old])
            for new, old in zip(expected_new_dec, gaia_ref_epoch_coords["dec"])
        ]
    )
    pmra = gaia_ref_epoch_coords["pmra"] / np.cos(np.deg2rad(average_dec))

    # angular displacement components
    gaia_ref_epoch_coords["pm_delta_dec"] = u.Quantity(
        (epoch - gaia_ref_epoch) * gaia_ref_epoch_coords["pmdec"] / 1000,
        unit=u.arcsec,
    ).to(u.deg)
    gaia_ref_epoch_coords["pm_delta_ra"] = u.Quantity(
        (epoch - gaia_ref_epoch) * (pmra / 1000), unit=u.arcsec
    ).to(u.deg)


def get_parallax_correction(epoch, gaia_ref_epoch_coords):
    """
    Calculates the parallax correction for a given epoch and Gaia reference epoch
    coordinates.

    Parameters
    ----------
    epoch : float
        The epoch for which to calculate the parallax correction.
    gaia_ref_epoch_coords : dict
        A dictionary containing the Gaia reference epoch coordinates:
        - "ra" : float
            The right ascension in degrees.
        - "dec" : float
            The declination in degrees.
        - "parallax" : float
            The parallax in milliarcseconds (mas).

    Returns
    -------
    None

    Notes
    -----
    This function calculates the parallax correction for a given epoch and Gaia
    reference epoch coordinates. It uses the `get_parallax_correction_barycenter`
    and `get_parallax_correction_mast` functions to obtain the parallax corrections
    based on different coordinate frames.

    Examples
    --------
    This function is typically used to add parallax correction columns to a main table
    of Gaia reference epoch coordinates.

    .. code-block:: python

        epoch = 2023.5
        gaia_coords = {
            "ra": 180.0,
            "dec": 30.0,
            "parallax": 2.5
        }
        get_parallax_correction(epoch, gaia_coords)
    """  # noqa: E501

    # get parallax correction using textbook calculations (i.e. Earth's barycenter)
    parallax_corr = get_parallax_correction_barycenter(
        epoch=epoch, gaia_ref_epoch_coords=gaia_ref_epoch_coords
    )

    # add parallax corrections columns to the main table
    gaia_ref_epoch_coords["parallax_delta_ra"] = parallax_corr[0]
    gaia_ref_epoch_coords["parallax_delta_dec"] = parallax_corr[1]


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
        l2.meta.target["proper_motion_epoch"] = "2016.0"
        # update wcsinfo
        update_wcsinfo(l2)
        # add a dummy WCS object
        create_wcs_for_tweakreg_pipeline(l2, shift_1=shift_1, shift_2=shift_2)
        l2_im = rdm.ImageModel(l2)
        return l2_im

    return _base_image


@pytest.mark.parametrize(
    "catalog, num_sources",
    [
        ("GAIADR1", 5),
        ("GAIADR2", 10),
        ("GAIADR3", 15),
    ],
)
def test_create_astrometric_catalog_variable_num_sources(
    tmp_path, catalog, num_sources, request
):
    """Test fetching data from supported catalogs with variable number of sources."""
    output_filename = "ref_cat.ecsv"
    img = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)
    res = create_astrometric_catalog(
        [img],
        catalog=catalog,
        output=os.path.join(tmp_path, output_filename),
        num_sources=num_sources,
    )

    assert len(res) == num_sources


def test_create_astrometric_catalog_write_results_to_disk(tmp_path, base_image):
    img = base_image(shift_1=1000, shift_2=1000)
    num_sources = 5
    output_filename = "output"

    # get list of supported write formats
    fh = StringIO()
    table.Table.write.list_formats(out=fh)
    fh.seek(0)
    list_of_supported_formats = [
        x.strip().split()[0]
        for x in fh.readlines()[2:]
        if x.strip().split()[1].lower() == "yes"
    ]
    # exclude data formats
    [
        list_of_supported_formats.remove(x)
        for x in [
            "asdf",
            "fits",
            "hdf5",
            "parquet",
            "pandas.html",
            "pandas.json",
            "pandas.csv",
        ]
    ]

    for table_format in list_of_supported_formats:
        res = create_astrometric_catalog(
            [img],
            catalog="GAIADR3",
            output=os.path.join(tmp_path, output_filename),
            table_format=table_format,
            num_sources=num_sources,
        )

        assert len(res) == num_sources
        assert os.path.exists(os.path.join(tmp_path, output_filename))


@pytest.mark.parametrize(
    "catalog, epoch",
    [
        ("GAIADR1", "2000.0"),
        ("GAIADR2", "2010"),
        ("GAIADR3", "2030.0"),
        ("GAIADR3", "J2000"),
        ("GAIADR3", 2030.0),
        ("GAIADR3", None),
    ],
)
def test_create_astrometric_catalog_using_epoch(tmp_path, catalog, epoch, request):
    """Test fetching data from supported catalogs for a specific epoch."""
    output_filename = "ref_cat.ecsv"
    img = request.getfixturevalue("base_image")(shift_1=1000, shift_2=1000)

    metadata_epoch = (
        epoch if epoch is not None else img.meta.target["proper_motion_epoch"]
    )
    metadata_epoch = float(
        "".join(c for c in str(metadata_epoch) if c == "." or c.isdigit())
    )

    res = create_astrometric_catalog(
        [img],
        catalog=catalog,
        output=os.path.join(tmp_path, output_filename),
        epoch=epoch,
    )

    assert np.equal(res["epoch"], float(metadata_epoch)).all()


def test_compute_radius():
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
    "ra, dec, sr, catalog_name",
    [
        (10, 10, 0.1, "GAIADR1"),
        (10, 10, 0.1, "GAIADR2"),
        (10, 10, 0.1, "GAIADR3"),
        (10, -10, 0.1, "GAIADR1"),
        (10, -10, 0.1, "GAIADR2"),
        (10, -10, 0.1, "GAIADR3"),
        (0, 0, 0.01, "GAIADR1"),
        (0, 0, 0.01, "GAIADR2"),
        (0, 0, 0.01, "GAIADR3"),
    ],
)
def test_get_catalog_using_valid_parameters(ra, dec, sr, catalog_name):
    """Test that get_catalog works properly with valid input parameters."""

    result = get_catalog(ra, dec, sr=sr, catalog=catalog_name)

    assert len(result) > 0


@pytest.mark.parametrize(
    "ra, dec, sr, catalog_name",
    [
        (10, 10, 0.1, "GAIDR3"),
        (-10, 10, 0.1, "GAIADR3"),
        (10, 100, 0.1, "GAIADR3"),
        (10, 100, 0.1, ""),
        (None, 100, 0.1, "GAIADR3"),
    ],
)
def test_get_catalog_using_invalid_parameters(ra, dec, sr, catalog_name):
    """Test that get_catalog catches all exceptions when invalid input is provided."""

    with pytest.raises(Exception) as exec_info:
        get_catalog(ra, dec, sr=sr, catalog=catalog_name)

    assert exec_info.typename.lower() == "exception"


def test_get_catalog_valid_parameters_but_no_sources_returned():
    """Test that get_catalog raises an exception if no sources are found."""

    with pytest.raises(Exception) as exec_info:
        # use a small search radius (0.5 arcsec) to force no source detection
        get_catalog(10, 10, sr=0.00014, catalog="GAIADR3")

    assert exec_info.typename.lower() == "exception"


@pytest.mark.parametrize(
    "ra, dec, epoch",
    [
        (10, 10, 2000),
        (10, 10, 2010.3),
        (10, 10, 2030),
        (10, -10, 2000),
        (10, -10, 2010.3),
        (10, -10, 2030),
        (0, 0, 2000),
        (0, 0, 2010.3),
        (0, 0, 2030),
    ],
)
def test_get_catalog_using_epoch(ra, dec, epoch):
    """Test that get_catalog returns coordinates corrected by proper motion
    and parallax. The idea is to fetch data for a specific epoch from the MAST VO API
    and compare them with the expected coordinates for that epoch.
    First, the data for a specific coordinates and epoch are fetched from the MAST VO
    API. Then, the data for the same coordinates are fetched for the Gaia's reference
    epoch of 2016.0, and corrected for proper motion and parallax using explicit
    calculations for the initially specified epoch. We then compare the results between
    the returned coordinates from the MAST VO API and the manually updated
    coordinates."""

    result = get_catalog(ra, dec, epoch=epoch)

    # updated coordinates at the provided epoch
    returned_ra = np.array(result["ra"])
    returned_dec = np.array(result["dec"])

    # get GAIA data and update coords to requested epoch using pm measurements
    gaia_ref_epoch = 2016.0
    gaia_ref_epoch_coords_all = get_catalog(ra, dec, epoch=gaia_ref_epoch)

    gaia_ref_epoch_coords = gaia_ref_epoch_coords_all  # [mask]

    # calculate proper motion corrections
    get_proper_motion_correction(
        epoch=epoch,
        gaia_ref_epoch_coords=gaia_ref_epoch_coords,
        gaia_ref_epoch=gaia_ref_epoch,
    )
    # calculate parallax corrections
    get_parallax_correction(epoch=epoch, gaia_ref_epoch_coords=gaia_ref_epoch_coords)

    # calculate the expected coordinates value after corrections have been applied to
    # Gaia's reference epoch coordinates

    # textbook (barycentric frame)
    expected_ra = (
        gaia_ref_epoch_coords["ra"]
        + gaia_ref_epoch_coords["pm_delta_ra"]
        + gaia_ref_epoch_coords["parallax_delta_ra"]
    )
    expected_dec = (
        gaia_ref_epoch_coords["dec"]
        + gaia_ref_epoch_coords["pm_delta_dec"]
        + gaia_ref_epoch_coords["parallax_delta_dec"]
    )

    assert len(result) > 0

    # adopted tolerance: 2.8e-9 deg -> 10 uas (~0.0001 pix)
    assert np.median(returned_ra - expected_ra) < 2.8e-9
    assert np.median(returned_dec - expected_dec) < 2.8e-9

    assert mad_std(returned_ra - expected_ra) < 2.8e-9
    assert mad_std(returned_dec - expected_dec) < 2.8e-9


def test_get_catalog_timeout():
    """Test that get_catalog can timeout."""

    with pytest.raises(Exception) as exec_info:
        for dt in np.arange(1, 0, -0.01):
            try:
                get_catalog(10, 10, sr=0.1, catalog="GAIADR3", timeout=dt)
            except requests.exceptions.ConnectionError:
                # ignore if it's a connection error instead of timeout
                pass

    assert exec_info.type == requests.exceptions.Timeout


def test_get_catalog_raises_connection_error(monkeypatch):
    """Test that get_catalog can raise a connection error."""

    monkeypatch.setattr("requests.get", MockConnectionError)

    with pytest.raises(Exception) as exec_info:
        get_catalog(10, 10, sr=0.1, catalog="GAIADR3")

    assert exec_info.type == requests.exceptions.ConnectionError
