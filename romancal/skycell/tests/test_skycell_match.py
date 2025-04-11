"""
Unit tests for proj_match.

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

import asdf
import astropy.coordinates as coord
import astropy.modeling.models as amm
import astropy.units as u
import gwcs.coordinate_frames as cf
import gwcs.wcs as wcs
import numpy as np
import pytest
import spherical_geometry.vector as sgv
from spherical_geometry.vector import rotate_around as rotate

import romancal.skycell.match as sm

SKYCELLS_TABLE_SUBSET_PATH = Path(__file__).parent / "skycells_subset.asdf"
# do not use load_patch_table here as it will modify global state
with asdf.open(SKYCELLS_TABLE_SUBSET_PATH) as skycells_subset:
    SKYCELLS_TABLE_SUBSET = skycells_subset.copy()

crecord = SKYCELLS_TABLE_SUBSET["roman"]["skycells"][0]
cra = crecord["ra_corn3"]
cdec = crecord["dec_corn3"]

cpa = 45.0
csize = 0.001
e = 0.0011  # epsilon offset in degrees


@pytest.fixture(autouse=True)
def override_patch_table(monkeypatch):
    """
    For the tests in this file, monkeypatch the global
    SKYCELLS_TABLE to a smaller SKYCELLS_TABLE_SUBSET to allow these tests
    to run without access to the full patch table.
    """
    monkeypatch.setattr(sm, "SKYCELLS_TABLE", SKYCELLS_TABLE_SUBSET)
    yield


def mk_im_corners(ra, dec, pa, size):
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
    radecrect = np.array(frect).transpose()
    return radecrect


def mk_gwcs(ra=cra, dec=cdec, pa=cpa, bounding_box=None, pixel_shape=None):
    """
    Construct a GWCS model for testing the patch matching when provided a WCS
    This just implements a basic tangent projection with specified ra, dec, and
    position angle
    """
    transform = (amm.Shift(-2048) & amm.Shift(-2048)) | (
        amm.Scale(0.11 / 3600.0) & amm.Scale(0.11 / 3600.0)
        | amm.Rotation2D(pa)
        | amm.Pix2Sky_TAN()
        | amm.RotateNative2Celestial(cra, cdec, 180.0)
    )
    detector_frame = cf.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = cf.CelestialFrame(
        reference_frame=coord.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = wcs.WCS([(detector_frame, transform), (sky_frame, None)])
    if pixel_shape is not None:
        wcsobj.pixel_shape = pixel_shape
    if bounding_box is not None:
        wcsobj.bounding_box = bounding_box
    return wcsobj


@pytest.mark.parametrize(
    "pars, expected",
    [
        (
            (cra, cdec + e, cpa, csize),
            (
                "225p90x26y50",
                "225p90x26y51",
            ),
        ),
        (
            (cra, cdec - e, cpa, csize),
            (
                "225p90x25y50",
                "225p90x25y51",
                "225p90x26y50",
                "225p90x26y51",
            ),
        ),
        (
            (cra + e, cdec, cpa, csize),
            (
                "225p90x25y50",
                "225p90x25y51",
                "225p90x26y50",
                "225p90x26y51",
            ),
        ),
        (
            (cra - e, cdec, cpa, csize),
            (
                "225p90x25y50",
                "225p90x25y51",
                "225p90x26y50",
                "225p90x26y51",
            ),
        ),
        (
            (cra, cdec, cpa, csize),
            (
                "225p90x25y50",
                "225p90x25y51",
                "225p90x26y50",
                "225p90x26y51",
            ),
        ),
        (
            (cra, cdec, cpa, 0.5),
            (
                "225p90x25y50",
                "225p90x25y51",
                "225p90x26y46",
                "225p90x26y47",
                "225p90x26y48",
                "225p90x26y49",
                "225p90x26y50",
                "225p90x26y51",
                "225p90x26y52",
                "225p90x26y53",
                "225p90x26y54",
                "225p90x26y55",
                "225p90x27y47",
                "225p90x27y48",
                "225p90x27y49",
                "225p90x27y50",
                "225p90x27y51",
                "225p90x27y52",
                "225p90x27y53",
                "225p90x27y54",
                "225p90x28y48",
                "225p90x28y49",
                "225p90x28y50",
                "225p90x28y51",
                "225p90x28y52",
                "225p90x28y53",
                "225p90x29y49",
                "225p90x29y50",
                "225p90x29y51",
                "225p90x29y52",
                "225p90x30y50",
                "225p90x30y51",
            ),
        ),
    ],
)
def test_corners(pars, expected):
    corners = mk_im_corners(*pars)
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(corners)
    # map matches to absolute index
    mmatches = tuple(
        [
            SKYCELLS_TABLE_SUBSET["roman"]["skycells"][index]["name"]
            for index in intersecting_skycells
        ]
    )
    assert tuple(mmatches) == expected


def test_wcs_corners():
    imshape = (4096, 4096)
    wcsobj = mk_gwcs()
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(
        wcsobj, image_shape=imshape
    )
    mmatches = tuple(
        [
            SKYCELLS_TABLE_SUBSET["roman"]["skycells"][index]["name"]
            for index in intersecting_skycells
        ]
    )
    assert tuple(mmatches) == (
        "225p90x25y50",
        "225p90x25y51",
        "225p90x26y49",
        "225p90x26y50",
        "225p90x26y51",
        "225p90x26y52",
        "225p90x27y50",
        "225p90x27y51",
    )

    wcsobj.pixel_shape = imshape
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
    mmatches = tuple(
        [
            SKYCELLS_TABLE_SUBSET["roman"]["skycells"][index]["name"]
            for index in intersecting_skycells
        ]
    )
    assert tuple(mmatches) == (
        "225p90x25y50",
        "225p90x25y51",
        "225p90x26y49",
        "225p90x26y50",
        "225p90x26y51",
        "225p90x26y52",
        "225p90x27y50",
        "225p90x27y51",
    )

    wcsobj.pixel_shape = None
    wcsobj.bounding_box = ((-0.5, 4096 - 0.5), (-0.5, 4096 - 0.5))
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
    mmatches = tuple(
        [
            SKYCELLS_TABLE_SUBSET["roman"]["skycells"][match]["name"]
            for match in intersecting_skycells
        ]
    )
    assert tuple(mmatches) == (
        "225p90x25y50",
        "225p90x25y51",
        "225p90x26y49",
        "225p90x26y50",
        "225p90x26y51",
        "225p90x26y52",
        "225p90x27y50",
        "225p90x27y51",
    )

    wcsobj.bounding_box = None
    with pytest.raises(ValueError):
        intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
