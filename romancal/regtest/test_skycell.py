"""Unit tests for sky cell WCS on global sky map"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from romancal.skycell import skymap
from romancal.skycell.tests.test_skycell import assert_corners_on_pixel_corners

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


def test_skycells():
    skycells = skymap.SkyCells.from_names(TEST_SKYCELLS)

    assert sorted(skycells.names) == sorted(TEST_SKYCELLS)

    assert skycells.radec_corners.shape == (len(TEST_SKYCELLS), 4, 2)
    assert skycells.vectorpoint_corners.shape == (len(TEST_SKYCELLS), 4, 3)

    assert skycells.radec_centers.shape == (len(TEST_SKYCELLS), 2)
    assert skycells.vectorpoint_centers.shape == (len(TEST_SKYCELLS), 3)

    assert len(skycells.polygons) == len(TEST_SKYCELLS)


def test_skycell_from_name():
    skycells = skymap.SkyCells.from_names(TEST_SKYCELLS)
    assert len(skycells) == len(TEST_SKYCELLS)

    with pytest.raises(KeyError):
        # this sky cell should not exist, even in the global skymap
        skymap.SkyCells.from_names(["270p65x99y70"])


@pytest.mark.parametrize("name", TEST_SKYCELLS)
def test_skycell_wcs_pixel_to_world(name):
    skycell = skymap.SkyCells.from_names([name])

    wcsobj = skycell.wcs[0]

    # the center pixel falls on the recorded center, whatever the handedness
    assert_allclose_lonlat(
        wcsobj(
            (skycell.pixel_shape[0] / 2.0) - 0.5,
            (skycell.pixel_shape[1] / 2.0) - 0.5,
        ),
        skycell.radec_centers[0],
        rtol=1e-7,
    )


@pytest.mark.parametrize("name", TEST_SKYCELLS)
def test_skycell_wcs_world_to_pixel(name):
    skycell = skymap.SkyCells.from_names([name])

    assert_corners_on_pixel_corners(
        skycell.wcs[0], skycell.radec_corners[0], skycell.pixel_shape
    )


def test_skycells_projection_regions():
    rng = np.random.default_rng()
    lon = rng.standard_normal(1)
    lon = lon / np.max(np.abs(lon)) * 180 + 180
    lat = rng.standard_normal(1)
    lat = lat / np.max(np.abs(lat)) * 90

    assert len(skymap.SKYMAP.skycells.projection_regions) == len(
        skymap.SKYMAP.model.skycells
    )


def test_skycells_containing_points():
    rng = np.random.default_rng()
    lon = rng.standard_normal(10000)
    lon = lon / np.max(np.abs(lon)) * 180 + 180
    lat = rng.standard_normal(10000)
    lat = lat / np.max(np.abs(lat)) * 90
    radec = np.stack([lon, lat], axis=1)

    skycells_containing_points = skymap.SKYMAP.skycells.containing(radec)
    point_indices_outside_skycells = [
        point_index
        for point_index in np.arange(radec.shape[0])
        if not any(
            point_index not in skycell_point_indices
            for skycell_point_indices in skycells_containing_points.values()
        )
    ]

    assert len(point_indices_outside_skycells) == 0, (
        f"{len(point_indices_outside_skycells)} / {radec.shape[0]} points do not lie within any skycell"
    )


def test_skycells_cores_containing_points():
    rng = np.random.default_rng()
    lon = rng.standard_normal(10000)
    lon = lon / np.max(np.abs(lon)) * 180 + 180
    lat = rng.standard_normal(10000)
    lat = lat / np.max(np.abs(lat)) * 90
    radec = np.stack([lon, lat], axis=1)

    skycells_exclusively_containing_points = skymap.SKYMAP.skycells.cores_containing(
        radec
    )
    point_indices_outside_core = [
        point_index
        for point_index in np.arange(radec.shape[0])
        if not any(
            point_index not in skycell_point_indices
            for skycell_point_indices in skycells_exclusively_containing_points.values()
        )
    ]

    # each point on the sphere MUST belong to exactly one skycell
    assert len(point_indices_outside_core) == 0, (
        f"{len(point_indices_outside_core)} / {radec.shape[0]} points do not lie within the exclusive zone of any skycell"
    )
