"""Unit tests for skycell WCS on global sky map"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import romancal.skycell.skymap as sc

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(
        np.array(actual) % 360, np.array(desired) % 360, rtol=rtol, atol=atol
    )


def test_skycell_from_name():
    with pytest.raises(KeyError):
        # this sky cell should not exist, even in the global sky map
        sc.SkyCell.from_name("270p65x99y70")


@pytest.mark.parametrize(
    "name",
    [
        "000p86x61y68",
        "045p86x29y34",
        # north pole
        "225p90x25y49",
        # south pole
        "135m90x40y46",
    ],
)
def test_skycell_wcs_pixel_to_world(name):
    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs

    # forward transform to radec corners
    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose_lonlat(
        np.array(
            wcs(
                *np.array(
                    [
                        (0.5, 0.5),
                        (skycell.pixel_shape[0] + 0.5, 0.5),
                        (skycell.pixel_shape[0] + 0.5, skycell.pixel_shape[1] + 0.5),
                        (0.5, skycell.pixel_shape[1] + 0.5),
                    ]
                ).T
            )
        ).T,
        skycell.radec_corners,
        rtol=1e-7,
    )


@pytest.mark.parametrize(
    "name",
    [
        "000p86x30y34",
        "045p86x29y34",
        # north pole
        "225p90x25y49",
        # south pole
        "135m90x40y46",
    ],
)
def test_skycell_wcs_world_to_pixel(name):
    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs

    # inverse transform to pixel corners
    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose(
        np.array(wcs.invert(*skycell.radec_corners.T, with_bounding_box=False)).T,
        (
            [
                (0.5, 0.5),
                (skycell.pixel_shape[0] + 0.5, 0.5),
                (skycell.pixel_shape[0] + 0.5, skycell.pixel_shape[1] + 0.5),
                (0.5, skycell.pixel_shape[1] + 0.5),
            ]
        ),
        rtol=1e-5,
    )
