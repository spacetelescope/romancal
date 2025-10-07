"""Unit tests for skycell WCS on global sky map"""

import numpy as np
import pytest
import spherical_geometry.vector as sgv
from numpy.testing import assert_allclose
from scipy.spatial import KDTree

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


def test_skycells():
    skycells = sc.SkyCells.from_names(TEST_SKYCELLS)

    assert sorted(skycells.names) == sorted(TEST_SKYCELLS)

    assert skycells.radec_corners.shape == (len(TEST_SKYCELLS), 4, 2)
    assert skycells.vectorpoint_corners.shape == (len(TEST_SKYCELLS), 4, 3)

    assert skycells.radec_centers.shape == (len(TEST_SKYCELLS), 2)
    assert skycells.vectorpoint_centers.shape == (len(TEST_SKYCELLS), 3)

    assert len(skycells.polygons) == len(TEST_SKYCELLS)


def test_skycells_core_contains_points():
    rng = np.random.default_rng()
    lon = rng.standard_normal(1000)
    lon = lon / np.max(np.abs(lon)) * 180 + 180
    lat = rng.standard_normal(1000)
    lat = lat / np.max(np.abs(lat)) * 90
    radec = np.stack([lon, lat], axis=1)

    skycells_containing_points = sc.SKYMAP.skycells_at(radec)
    point_indices_outside_skycells = [
        point_index
        for point_index, skycell_indices in enumerate(skycells_containing_points)
        if len(skycell_indices) == 0
    ]

    assert len(point_indices_outside_skycells) == 0, (
        f"{len(point_indices_outside_skycells)} / {radec.shape[0]} points do not lie within any skycell"
    )

    point_indices_outside_core = [
        point_index
        for point_index, skycells in enumerate(skycells_containing_points)
        if not len(
            np.where(
                [skycell.core_contains(radec[point_index, :]) for skycell in skycells]
            )[0]
        )
        == 1
    ]

    # each point on the sphere MUST belong to exactly one skycell
    assert len(point_indices_outside_core) == 0, (
        f"{len(point_indices_outside_core)} / {radec.shape[0]} points do not lie within the exclusive zone of any skycell"
    )
