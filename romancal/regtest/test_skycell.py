"""Unit tests for skycell WCS on global sky map"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import romancal.skycell.skymap as sc

# mark all tests in this module
pytestmark = [pytest.mark.bigdata]

TEST_SKYCELLS = [
    "000p86x68y61",
    "045p86x34y29",
    # north pole
    "135p90x49y25",
    # south pole
    "225m90x46y40",
]


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(
        np.array(actual) % 360, np.array(desired) % 360, rtol=rtol, atol=atol
    )


def test_skycell_from_name():
    with pytest.raises(KeyError):
        # this sky cell should not exist, even in the global sky map
        sc.SkyCell.from_name("270p65x99y70")


@pytest.mark.parametrize("name", TEST_SKYCELLS)
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
                        (-0.5, -0.5),
                        (skycell.pixel_shape[0] - 0.5, -0.5),
                        (skycell.pixel_shape[0] - 0.5, skycell.pixel_shape[1] - 0.5),
                        (-0.5, skycell.pixel_shape[1] - 0.5),
                    ]
                ).T,
                with_bounding_box=False,
            )
        ).T,
        skycell.radec_corners,
        rtol=1e-7,
    )


@pytest.mark.parametrize("name", TEST_SKYCELLS)
def test_skycell_wcs_world_to_pixel(name):
    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs

    # inverse transform to pixel corners
    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose(
        np.array(wcs.invert(*skycell.radec_corners.T, with_bounding_box=False)).T,
        (
            [
                (-0.5, -0.5),
                (skycell.pixel_shape[0] - 0.5, -0.5),
                (skycell.pixel_shape[0] - 0.5, skycell.pixel_shape[1] - 0.5),
                (-0.5, skycell.pixel_shape[1] - 0.5),
            ]
        ),
        rtol=1e-5,
    )


def test_skycells(skymap_subset):
    skycells = sc.SkyCells.from_names(TEST_SKYCELLS, skymap=skymap_subset)

    assert sorted(skycells.names) == sorted(TEST_SKYCELLS)

    assert skycells.radec_corners.shape == (len(TEST_SKYCELLS), 4, 2)
    assert skycells.vectorpoint_corners.shape == (len(TEST_SKYCELLS), 4, 3)

    assert skycells.radec_centers.shape == (len(TEST_SKYCELLS), 2)
    assert skycells.vectorpoint_centers.shape == (len(TEST_SKYCELLS), 3)

    assert len(skycells.polygons)
