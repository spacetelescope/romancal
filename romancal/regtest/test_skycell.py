"""Unit tests for sky cell WCS on global sky map"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from romancal.skycell import skymap

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

    # forward transform to RA Dec corners
    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose_lonlat(
        np.array(
            wcsobj(
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
        skycell.radec_corners[0],
        rtol=1e-7,
    )


@pytest.mark.parametrize("name", TEST_SKYCELLS)
def test_skycell_wcs_world_to_pixel(name):
    skycell = skymap.SkyCells.from_names([name])

    wcsobj = skycell.wcs[0]

    # inverse transform to pixel corners
    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose(
        np.array(wcsobj.invert(*skycell.radec_corners.T, with_bounding_box=False)).T,
        [
            [
                (-0.5, -0.5),
                (skycell.pixel_shape[0] - 0.5, -0.5),
                (skycell.pixel_shape[0] - 0.5, skycell.pixel_shape[1] - 0.5),
                (-0.5, skycell.pixel_shape[1] - 0.5),
            ]
        ],
        rtol=1e-5,
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
    npoints = 100000
    radec = np.stack(
        [
            rng.uniform(0, 360, npoints),
            np.degrees(np.arcsin(rng.uniform(-1, 1, npoints))),
        ],
        axis=1,
    )

    skycells_containing_points = skymap.SKYMAP.skycells.containing(radec)
    skycells_exclusively_containing_points = skymap.SKYMAP.skycells.cores_containing(
        radec
    )
    containing_pairs = {
        (skycell_index, point_index)
        for skycell_index, point_indices in skycells_containing_points.items()
        for point_index in point_indices
    }
    skycells_per_point = np.bincount(
        [point_index for _, point_index in containing_pairs], minlength=npoints
    )

    assert np.all(skycells_per_point >= 1), (
        f"{np.sum(skycells_per_point == 0)} / {npoints} points do not lie within any skycell"
    )
    # the skycell whose core contains a point also contains it
    assert all(
        (skycell_index, point_index) in containing_pairs
        for skycell_index, point_indices in skycells_exclusively_containing_points.items()
        for point_index in point_indices
    )


def test_skycells_cores_containing_points():
    rng = np.random.default_rng()
    npoints = 100000
    radec = np.stack(
        [
            rng.uniform(0, 360, npoints),
            np.degrees(np.arcsin(rng.uniform(-1, 1, npoints))),
        ],
        axis=1,
    )

    skycells_exclusively_containing_points = skymap.SKYMAP.skycells.cores_containing(
        radec
    )
    skycells_per_point = np.bincount(
        np.concatenate(list(skycells_exclusively_containing_points.values())),
        minlength=npoints,
    )

    # each point on the sphere MUST belong to exactly one skycell
    assert np.all(skycells_per_point == 1), (
        f"{np.sum(skycells_per_point == 0)} / {npoints} points lie within no skycell core, "
        f"{np.sum(skycells_per_point > 1)} within more than one"
    )
