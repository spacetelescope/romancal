"""
Unit tests for patch_match.

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

import inspect
import os
import os.path

import astropy.coordinates as coord
import astropy.modeling.models as amm
import astropy.units as u
import gwcs.coordinate_frames as cf
import gwcs.wcs as wcs
import numpy as np
import pytest
import spherical_geometry.vector as sgv

import romancal.patch_match.patch_match as pm

mpath = inspect.stack()[0][1]
pm.load_patch_table(os.path.join(os.path.split(mpath)[0], "patches_subset.asdf"))
rotate = sgv.rotate_around
absindex = 925050
patchtable = pm.PATCH_TABLE
crecord = patchtable[np.where(patchtable[:]["index"] == absindex)]
cra = crecord["ra_corn3"]
if len(cra) == 1:
    cra = cra[0]
cdec = crecord["dec_corn3"]
if len(cdec) == 1:
    cdec = cdec[0]
cpa = 45.0
csize = 0.001
e = 0.0011  # epsilon offset in degrees


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
        ((cra, cdec + e, cpa, csize), (925051, 925151)),
        ((cra, cdec - e, cpa, csize), (925050, 925150)),
        ((cra + e, cdec, cpa, csize), (925150, 925151)),
        ((cra - e, cdec, cpa, csize), (925050, 925051)),
        ((cra, cdec, cpa, csize), (925050, 925051, 925150, 925151)),
        (
            (cra, cdec, cpa, 0.5),
            (
                924750,
                924751,
                924849,
                924850,
                924851,
                924852,
                924948,
                924949,
                924950,
                924951,
                924952,
                924953,
                925047,
                925048,
                925049,
                925050,
                925051,
                925052,
                925053,
                925054,
                925147,
                925148,
                925149,
                925150,
                925151,
                925152,
                925153,
                925154,
                925248,
                925249,
                925250,
                925251,
                925252,
                925253,
                925349,
                925350,
                925351,
                925352,
                925450,
                925451,
            ),
        ),
    ],
)
def test_corners(pars, expected):
    corners = mk_im_corners(*pars)
    matches, close = pm.find_patch_matches(corners)
    # map matches to absolute index
    mmatches = tuple([patchtable[match]["index"] for match in matches])
    assert tuple(mmatches) == expected


def test_wcs_corners():
    imshape = (4096, 4096)
    wcsobj = mk_gwcs()
    matches, close = pm.find_patch_matches(wcsobj, image_shape=imshape)
    mmatches = tuple([patchtable[match]["index"] for match in matches])
    assert tuple(mmatches) == (925050, 925051, 925150, 925151)
    wcsobj.pixel_shape = imshape
    matches, close = pm.find_patch_matches(wcsobj)
    mmatches = tuple([patchtable[match]["index"] for match in matches])
    assert tuple(mmatches) == (925050, 925051, 925150, 925151)
    wcsobj.pixel_shape = None
    wcsobj.bounding_box = ((-0.5, 4096 - 0.5), (-0.5, 4096 - 0.5))
    matches, close = pm.find_patch_matches(wcsobj)
    mmatches = tuple([patchtable[match]["index"] for match in matches])
    assert tuple(mmatches) == (925050, 925051, 925150, 925151)
    wcsobj.bounding_box = None
    with pytest.raises(ValueError):
        matches, close = pm.find_patch_matches(wcsobj)
