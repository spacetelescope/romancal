"""
Unit tests for skycell.match.

These tests depend very strongly on the contents of the referenced table of patches.
Changes to the contents of this table will require changes to the tests for
any unit tests that depend on specific matches to patches in the table.

Any changes to the matching algorithm should be completely separate from any
changes to the table contents.

The tests include the following cases to validate the matches
1) Simple case at ra, dec, pa = 0, 0, 0
2) 1) translated to ra, de, pa = 180, 40, 0
3) Same as 2) but with pa = 45 and 60.
4) Same as 3) (pa=45 case) but with the lower corner just above, below, to the right
   and to the left of a 4 patch corner
   (assuming non-overlapping patchs within a common tangent point).
   This requires identifying the ra, dec of such a corner in the table.
5) A test of a WCS provided example.

Most of these tests check to see if the matches are what are expected, by index of
the table (the table format is expected to include the index of the entry as one
of its columns so that subsets of the table selected to reduce the filesize still
retain the same index obtained.)
"""

from pathlib import Path

import astropy.coordinates as coord
import astropy.modeling.models as amm
import astropy.units as u
import numpy as np
import pytest
import spherical_geometry.vector as sgv
from gwcs import WCS, coordinate_frames

import romancal.skycell.match as sm
import romancal.skycell.skymap as sc

TEST_POINTS = [
    (0.88955854, 87.53857137),
    (20.6543883, 87.60498618),
    (343.19474696, 85.05565535),
    (8.94286202, 85.50465173),
    (27.38417684, 85.03404907),
    (310.53503934, 88.56749324),
]
EPSILON = 0.0011  # epsilon offset in degrees
DATA_DIRECTORY = Path(__file__).parent / "data"


@pytest.fixture(autouse=True)
def override_skymap(monkeypatch):
    """
    For the tests in this file, monkeypatch the global
    skymap path to a smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    monkeypatch.setattr(sc.SKYMAP, "path", DATA_DIRECTORY / "skymap_subset.asdf")
    yield


def mk_im_corners(
    ra: float, dec: float, pa: float, size: float
) -> list[tuple[float, float]]:
    """
    Generate 4 image corners of a square with the center at the supplied
    side size, ra, dec, and position angle (all in degrees).
    """
    # Generate 4 unit vectors at ra, dec = (0 , 0)
    center = sgv.lonlat_to_vector(0.0, 0.0)
    radecvec = sgv.lonlat_to_vector(ra, dec)
    zaxis = (0.0, 0.0, 1.0)
    yaxis = (0.0, 1.0, 0.0)
    pp = sgv.rotate_around(
        *(sgv.rotate_around(*(center + yaxis + (-size / 2,))) + zaxis + (+size / 2,))
    )
    pm = sgv.rotate_around(
        *(sgv.rotate_around(*(center + yaxis + (+size / 2,))) + zaxis + (+size / 2,))
    )
    mp = sgv.rotate_around(
        *(sgv.rotate_around(*(center + yaxis + (-size / 2,))) + zaxis + (-size / 2,))
    )
    mm = sgv.rotate_around(
        *(sgv.rotate_around(*(center + yaxis + (+size / 2,))) + zaxis + (-size / 2,))
    )
    rect = [pp, mp, mm, pm]

    # Now move to requested ra and dec
    trect = [
        sgv.rotate_around(
            *(sgv.rotate_around(*(vec + yaxis + (-dec,))) + zaxis + (ra,))
        )
        for vec in rect
    ]
    # Rotate to desired position angle
    rrect = [sgv.rotate_around(*(vec + radecvec + (pa,))) for vec in trect]
    frect = [sgv.vector_to_lonlat(*vec) for vec in rrect]
    # Reorganize by ra, dec arrays
    radecrect = np.array(frect)
    return radecrect


def mk_gwcs(ra, dec, pa, bounding_box=None, pixel_shape=None) -> WCS:
    """
    Construct a GWCS model for testing the patch matching when provided a WCS
    This just implements a basic tangent projection with specified ra, dec, and
    position angle
    """
    transform = (amm.Shift(-2048) & amm.Shift(-2048)) | (
        amm.Scale(0.11 / 3600.0) & amm.Scale(0.11 / 3600.0)
        | amm.Rotation2D(pa)
        | amm.Pix2Sky_TAN()
        | amm.RotateNative2Celestial(ra, dec, 180.0)
    )
    detector_frame = coordinate_frames.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = coordinate_frames.CelestialFrame(
        reference_frame=coord.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = WCS([(detector_frame, transform), (sky_frame, None)])
    if pixel_shape is not None:
        wcsobj.pixel_shape = pixel_shape
    if bounding_box is not None:
        wcsobj.bounding_box = bounding_box
    return wcsobj


@pytest.mark.parametrize(
    "test_point,offset,rotation,size,expected_skycell_names",
    [
        (
            TEST_POINTS[0],
            (0, 0),
            45,
            0.001,
            (
                "000p86x50y65",
                "000p86x50y66",
                "000p86x51y65",
                "000p86x51y66",
            ),
        ),
        (
            TEST_POINTS[0],
            (0, +EPSILON),
            45,
            0.001,
            (
                "000p86x50y66",
                "000p86x51y66",
            ),
        ),
        (
            TEST_POINTS[1],
            (0, +EPSILON),
            45,
            0.001,
            (
                "000p86x61y69",
                "000p86x62y69",
            ),
        ),
        (
            TEST_POINTS[1],
            (0, -EPSILON),
            45,
            0.001,
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            TEST_POINTS[1],
            (+EPSILON, 0),
            45,
            0.001,
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            TEST_POINTS[1],
            (-EPSILON, 0),
            45,
            0.001,
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            TEST_POINTS[1],
            (0, 0),
            45,
            0.001,
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            TEST_POINTS[0],
            (0, 0),
            45,
            0.3,
            (
                "000p86x48y65",
                "000p86x48y66",
                "000p86x49y64",
                "000p86x49y65",
                "000p86x49y66",
                "000p86x49y67",
                "000p86x50y63",
                "000p86x50y64",
                "000p86x50y65",
                "000p86x50y66",
                "000p86x50y67",
                "000p86x50y68",
                "000p86x51y63",
                "000p86x51y64",
                "000p86x51y65",
                "000p86x51y66",
                "000p86x51y67",
                "000p86x51y68",
                "000p86x52y64",
                "000p86x52y65",
                "000p86x52y66",
                "000p86x52y67",
                "000p86x53y65",
                "000p86x53y66",
            ),
        ),
        (
            TEST_POINTS[1],
            (0, 0),
            45,
            0.5,
            (
                "000p86x57y67",
                "000p86x57y68",
                "000p86x58y66",
                "000p86x58y67",
                "000p86x58y68",
                "000p86x58y69",
                "000p86x58y70",
                "000p86x59y66",
                "000p86x59y67",
                "000p86x59y68",
                "000p86x59y69",
                "000p86x59y70",
                "000p86x59y71",
                "000p86x59y72",
                "000p86x60y65",
                "000p86x60y66",
                "000p86x60y67",
                "000p86x60y68",
                "000p86x60y69",
                "000p86x60y70",
                "000p86x60y71",
                "000p86x60y72",
                "000p86x60y73",
                "000p86x61y65",
                "000p86x61y66",
                "000p86x61y67",
                "000p86x61y68",
                "000p86x61y69",
                "000p86x61y70",
                "000p86x61y71",
                "000p86x61y72",
                "000p86x61y73",
                "000p86x62y64",
                "000p86x62y65",
                "000p86x62y66",
                "000p86x62y67",
                "000p86x62y68",
                "000p86x62y69",
                "000p86x62y70",
                "000p86x62y71",
                "000p86x63y64",
                "000p86x63y65",
                "000p86x63y66",
                "000p86x63y67",
                "000p86x63y68",
                "000p86x63y69",
                "000p86x64y65",
                "000p86x64y66",
                "045p86x36y65",
                "045p86x36y66",
                "045p86x37y66",
                "045p86x37y67",
                "045p86x37y68",
                "045p86x37y69",
                "045p86x38y66",
                "045p86x38y67",
                "045p86x38y68",
                "045p86x38y69",
                "045p86x38y70",
                "045p86x38y71",
                "045p86x39y66",
                "045p86x39y67",
                "045p86x39y68",
                "045p86x39y69",
                "045p86x39y70",
                "045p86x39y71",
                "045p86x39y72",
                "045p86x39y73",
                "045p86x39y74",
                "045p86x40y67",
                "045p86x40y68",
                "045p86x40y69",
                "045p86x40y70",
                "045p86x40y71",
                "045p86x41y67",
                "045p86x41y68",
                "045p86x41y69",
            ),
        ),
        (
            TEST_POINTS[2],
            (0, 0),
            0,
            0.4,
            (
                "000p86x27y33",
                "000p86x27y34",
                "000p86x28y32",
                "000p86x28y33",
                "000p86x28y34",
                "000p86x28y35",
                "000p86x28y36",
                "000p86x28y37",
                "000p86x29y32",
                "000p86x29y33",
                "000p86x29y34",
                "000p86x29y35",
                "000p86x29y36",
                "000p86x29y37",
                "000p86x29y38",
                "000p86x30y32",
                "000p86x30y33",
                "000p86x30y34",
                "000p86x30y35",
                "000p86x30y36",
                "000p86x30y37",
                "000p86x30y38",
                "000p86x31y31",
                "000p86x31y32",
                "000p86x31y33",
                "000p86x31y34",
                "000p86x31y35",
                "000p86x31y36",
                "000p86x31y37",
                "000p86x32y31",
                "000p86x32y32",
                "000p86x32y33",
                "000p86x32y34",
                "000p86x32y35",
                "000p86x32y36",
                "000p86x32y37",
                "000p86x33y32",
                "000p86x33y33",
                "000p86x33y34",
                "000p86x33y35",
                "000p86x33y36",
                "000p86x33y37",
                "000p86x34y35",
                "000p86x34y36",
                "000p86x34y37",
            ),
        ),
        (
            TEST_POINTS[3],
            (-0.5, -0.5),
            0,
            0.001,
            ("000p86x60y32",),
        ),
        (
            TEST_POINTS[4],
            (0, 0),
            -62,
            0.2,
            (
                "045p86x28y33",
                "045p86x28y34",
                "045p86x28y35",
                "045p86x28y36",
                "045p86x29y33",
                "045p86x29y34",
                "045p86x29y35",
                "045p86x29y36",
                "045p86x30y33",
                "045p86x30y34",
                "045p86x30y35",
                "045p86x30y36",
                "045p86x31y33",
                "045p86x31y34",
                "045p86x31y35",
                "045p86x31y36",
            ),
        ),
        (
            TEST_POINTS[5],
            (0, 0),
            188,
            0.25,
            (
                "225p90x29y50",
                "225p90x29y51",
                "225p90x29y52",
                "225p90x29y53",
                "225p90x30y50",
                "225p90x30y51",
                "225p90x30y52",
                "225p90x30y53",
                "225p90x31y50",
                "225p90x31y51",
                "225p90x31y52",
                "225p90x31y53",
                "225p90x32y50",
                "225p90x32y51",
                "225p90x32y52",
                "225p90x32y53",
            ),
        ),
    ],
)
def test_skycell_match(test_point, offset, rotation, size, expected_skycell_names):
    corners = mk_im_corners(*test_point + np.array(offset), rotation, size)

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(corners)

    skycell_names = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()

    assert sorted(skycell_names) == sorted(expected_skycell_names)


@pytest.mark.parametrize(
    "test_point,expected_skycell_names",
    [
        (
            TEST_POINTS[1],
            [
                "000p86x62y69",
                "000p86x62y68",
                "000p86x61y69",
                "000p86x61y68",
                "000p86x63y69",
                "000p86x61y70",
                "000p86x62y67",
                "000p86x60y68",
                "045p86x37y69",
                "045p86x37y68",
                "045p86x38y69",
            ],
        )
    ],
)
def test_match_from_wcs(test_point, expected_skycell_names):
    wcsobj = mk_gwcs(*test_point, 45)
    imshape = (4096, 4096)

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(
        wcsobj, image_shape=imshape
    )

    skycell_names = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()

    assert skycell_names == expected_skycell_names


@pytest.mark.parametrize(
    "test_point,expected_skycell_names",
    [
        (
            TEST_POINTS[1],
            [
                "000p86x62y69",
                "000p86x62y68",
                "000p86x61y69",
                "000p86x61y68",
                "000p86x63y69",
                "000p86x61y70",
                "000p86x62y67",
                "000p86x60y68",
                "045p86x37y69",
                "045p86x37y68",
                "045p86x38y69",
            ],
        )
    ],
)
def test_match_from_wcs_with_imshape(test_point, expected_skycell_names):
    wcsobj = mk_gwcs(*test_point, 45)
    imshape = (4096, 4096)
    wcsobj.pixel_shape = imshape

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)

    skycell_names = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()

    assert skycell_names == expected_skycell_names


@pytest.mark.parametrize(
    "test_point,expected_skycell_names",
    [
        (
            TEST_POINTS[1],
            [
                "000p86x62y69",
                "000p86x62y68",
                "000p86x61y69",
                "000p86x61y68",
                "000p86x63y69",
                "000p86x61y70",
                "000p86x62y67",
                "000p86x60y68",
                "045p86x37y69",
                "045p86x37y68",
                "045p86x38y69",
            ],
        )
    ],
)
def test_match_from_wcs_with_bbox(test_point, expected_skycell_names):
    wcsobj = mk_gwcs(*test_point, 45)
    wcsobj.bounding_box = ((-0.5, 4096 - 0.5), (-0.5, 4096 - 0.5))

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)

    skycell_names = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()

    assert skycell_names == expected_skycell_names


@pytest.mark.parametrize("test_point", [TEST_POINTS[1]])
def test_match_from_wcs_without_imshape_or_bbox(test_point):
    wcsobj = mk_gwcs(*test_point, 45)
    wcsobj.bounding_box = None

    with pytest.raises(ValueError):
        intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
