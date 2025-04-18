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

import os
from pathlib import Path

import astropy.coordinates as coord
import astropy.modeling.models as amm
import astropy.units as u
import gwcs.coordinate_frames as cf
import gwcs.wcs as wcs
import numpy as np
import pytest
import spherical_geometry.vector as sgv
from spherical_geometry.vector import rotate_around as rotate

os.environ["SKYMAP_PATH"] = str(Path(__file__).parent / "skymap_subset.asdf")

import romancal.skycell.match as sm
import romancal.skycell.skymap as sc

crecord = sc.SKYMAP.skycells[0]
cra = crecord["ra_corn3"]
cdec = crecord["dec_corn3"]

cpa = 45.0
e = 0.0011  # epsilon offset in degrees


def mk_im_corners(
    ra: float, dec: float, pa: float, size: float
) -> tuple[
    tuple[float, float], tuple[float, float], tuple[float, float], tuple[float, float]
]:
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
            (cra, cdec + e, cpa, 0.001),
            (
                "315p86x50y75",
                "315p86x51y75",
            ),
        ),
        (
            (cra, cdec - e, cpa, 0.001),
            (
                "315p86x50y75",
                "315p86x51y75",
            ),
        ),
        (
            (cra + e, cdec, cpa, 0.001),
            (
                "315p86x50y75",
                "315p86x51y75",
            ),
        ),
        (
            (cra - e, cdec, cpa, 0.001),
            (
                "315p86x50y75",
                "315p86x51y75",
            ),
        ),
        (
            (cra, cdec, cpa, 0.001),
            (
                "315p86x50y75",
                "315p86x51y75",
            ),
        ),
        (
            (cra, cdec, cpa, 0.5),
            (
                "315p86x46y74",
                "315p86x46y75",
                "315p86x47y73",
                "315p86x47y74",
                "315p86x47y75",
                "315p86x48y72",
                "315p86x48y73",
                "315p86x48y74",
                "315p86x48y75",
                "315p86x49y71",
                "315p86x49y72",
                "315p86x49y73",
                "315p86x49y74",
                "315p86x49y75",
                "315p86x50y70",
                "315p86x50y71",
                "315p86x50y72",
                "315p86x50y73",
                "315p86x50y74",
                "315p86x50y75",
                "315p86x51y70",
                "315p86x51y71",
                "315p86x51y72",
                "315p86x51y73",
                "315p86x51y74",
                "315p86x51y75",
                "315p86x52y71",
                "315p86x52y72",
                "315p86x52y73",
                "315p86x52y74",
                "315p86x52y75",
                "315p86x53y72",
                "315p86x53y73",
                "315p86x53y74",
                "315p86x53y75",
                "315p86x54y73",
                "315p86x54y74",
                "315p86x54y75",
                "315p86x55y74",
                "315p86x55y75",
            ),
        ),
    ],
)
def test_skycell_match(pars, expected):
    corners = mk_im_corners(*pars)
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(corners)
    # map matches to absolute index
    mmatches = tuple(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    )
    assert mmatches == expected


def test_wcs_corners():
    imshape = (4096, 4096)
    wcsobj = mk_gwcs()
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(
        wcsobj, image_shape=imshape
    )
    mmatches = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()
    assert mmatches == [
        "315p86x49y74",
        "315p86x49y75",
        "315p86x50y73",
        "315p86x50y74",
        "315p86x50y75",
        "315p86x51y73",
        "315p86x51y74",
        "315p86x51y75",
        "315p86x52y74",
        "315p86x52y75",
    ]

    wcsobj.pixel_shape = imshape
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
    mmatches = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()
    assert mmatches == [
        "315p86x49y74",
        "315p86x49y75",
        "315p86x50y73",
        "315p86x50y74",
        "315p86x50y75",
        "315p86x51y73",
        "315p86x51y74",
        "315p86x51y75",
        "315p86x52y74",
        "315p86x52y75",
    ]

    wcsobj.pixel_shape = None
    wcsobj.bounding_box = ((-0.5, 4096 - 0.5), (-0.5, 4096 - 0.5))
    intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
    mmatches = np.array(
        [sc.SKYMAP.skycells[index]["name"] for index in intersecting_skycells]
    ).tolist()
    assert mmatches == [
        "315p86x49y74",
        "315p86x49y75",
        "315p86x50y73",
        "315p86x50y74",
        "315p86x50y75",
        "315p86x51y73",
        "315p86x51y74",
        "315p86x51y75",
        "315p86x52y74",
        "315p86x52y75",
    ]

    wcsobj.bounding_box = None
    with pytest.raises(ValueError):
        intersecting_skycells, nearby_skycells = sm.find_skycell_matches(wcsobj)
