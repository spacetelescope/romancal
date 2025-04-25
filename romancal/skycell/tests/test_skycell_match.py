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
from spherical_geometry.vector import rotate_around as rotate

import romancal.skycell.match as sm
import romancal.skycell.skymap as sc

e = 0.0011  # epsilon offset in degrees
cpa = 45.0


@pytest.fixture(autouse=True)
def override_skymap(monkeypatch):
    """
    For the tests in this file, monkeypatch the global
    skymap path to a smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    monkeypatch.setattr(sc.SKYMAP, "path", Path(__file__).parent / "skymap_subset.asdf")
    yield


@pytest.fixture
def sample_point() -> list[float]:
    """sample point of skycell 000p86x61y68, near the boundary of projection regions 1 and 2"""
    crecord = sc.SKYMAP.skycells[3560]
    return [crecord["ra_corn3"], crecord["dec_corn3"]]


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
    pp = rotate(*(rotate(*(center + yaxis + (-size / 2,))) + zaxis + (+size / 2,)))
    pm = rotate(*(rotate(*(center + yaxis + (+size / 2,))) + zaxis + (+size / 2,)))
    mp = rotate(*(rotate(*(center + yaxis + (-size / 2,))) + zaxis + (-size / 2,)))
    mm = rotate(*(rotate(*(center + yaxis + (+size / 2,))) + zaxis + (-size / 2,)))
    rect = [pp, mp, mm, pm]

    # Now move to requested ra and dec
    trect = [
        rotate(*(rotate(*(vec + yaxis + (-dec,))) + zaxis + (ra,))) for vec in rect
    ]
    # Rotate to desired position angle
    rrect = [rotate(*(vec + radecvec + (pa,))) for vec in trect]
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


@pytest.fixture
def sample_wcs(sample_point) -> WCS:
    return mk_gwcs(sample_point[0], sample_point[1], cpa)


@pytest.mark.parametrize(
    "pars, expected_skycell_names",
    [
        (
            (0, +e, cpa, 0.001),
            (
                "000p86x61y69",
                "000p86x62y69",
            ),
        ),
        (
            (0, -e, cpa, 0.001),
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            (+e, 0, cpa, 0.001),
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            (-e, 0, cpa, 0.001),
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            (0, 0, cpa, 0.001),
            (
                "000p86x61y68",
                "000p86x61y69",
                "000p86x62y68",
                "000p86x62y69",
            ),
        ),
        (
            (0, 0, cpa, 0.5),
            (
                "000p86x60y69",
                "000p86x61y68",
                "000p86x61y69",
                "000p86x61y70",
                "000p86x62y67",
                "000p86x62y68",
                "000p86x62y69",
                "000p86x62y70",
                "000p86x63y68",
                "000p86x63y69",
                "045p86x37y67",
                "045p86x37y68",
                "045p86x37y69",
                "045p86x38y67",
                "045p86x38y68",
                "045p86x38y69",
                "045p86x38y70",
                "045p86x38y71",
                "045p86x39y69",
                "045p86x39y70",
            ),
        ),
    ],
)
def test_skycell_match(pars, expected_skycell_names, sample_point):
    sample_point[0] += pars[0]
    sample_point[1] += pars[1]
    corners = mk_im_corners(*sample_point, pars[2], pars[3])
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(corners)
    # map matches to absolute index
    skycell_names = tuple(
        np.array(
            [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
        ).tolist()
    )
    assert sorted(skycell_names) == sorted(expected_skycell_names)


def test_match_from_wcs(sample_wcs):
    imshape = (4096, 4096)

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(
        sample_wcs, image_shape=imshape
    )
    skycell_names = [
        str(sc.SKYMAP.skycells[index]["name"]) for index in intersecting_skycells
    ]

    assert skycell_names == [
        "000p86x62y69",
        "000p86x62y68",
        "000p86x61y69",
        "000p86x61y68",
        "000p86x63y69",
        "000p86x61y70",
        "000p86x62y67",
        "045p86x37y69",
        "045p86x37y68",
        "045p86x38y69",
    ]


def test_match_from_wcs_with_imshape(sample_wcs):
    imshape = (4096, 4096)
    sample_wcs.pixel_shape = imshape

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(sample_wcs)
    skycell_names = [
        str(sc.SKYMAP.skycells[index]["name"]) for index in intersecting_skycells
    ]

    assert skycell_names == [
        "000p86x62y69",
        "000p86x62y68",
        "000p86x61y69",
        "000p86x61y68",
        "000p86x63y69",
        "000p86x61y70",
        "000p86x62y67",
        "045p86x37y69",
        "045p86x37y68",
        "045p86x38y69",
    ]


def test_match_from_wcs_with_bbox(sample_wcs):
    sample_wcs.bounding_box = ((-0.5, 4096 - 0.5), (-0.5, 4096 - 0.5))

    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(sample_wcs)
    skycell_names = [
        str(sc.SKYMAP.skycells[index]["name"]) for index in intersecting_skycells
    ]
    assert skycell_names == [
        "000p86x62y69",
        "000p86x62y68",
        "000p86x61y69",
        "000p86x61y68",
        "000p86x63y69",
        "000p86x61y70",
        "000p86x62y67",
        "045p86x37y69",
        "045p86x37y68",
        "045p86x38y69",
    ]


def test_match_from_wcs_without_imshape_or_bbox(sample_wcs):
    sample_wcs.bounding_box = None
    with pytest.raises(ValueError):
        intersecting_skycells, nearby_skycells = sm.find_skycell_matches(sample_wcs)
