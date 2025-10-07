"""Unit tests for skycell WCS on global sky map"""

import numpy as np
import pytest
import spherical_geometry.vector as sgv
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
        "000p86x68y61",
        "045p86x34y29",
        # north pole
        "135p90x49y25",
        # south pole
        "225m90x46y40",
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


@pytest.mark.parametrize(
    "name",
    [
        "000p86x34y30",
        "045p86x34y29",
        # north pole
        "135p90x49y25",
        # south pole
        "225m90x46y40",
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
                (-0.5, -0.5),
                (skycell.pixel_shape[0] - 0.5, -0.5),
                (skycell.pixel_shape[0] - 0.5, skycell.pixel_shape[1] - 0.5),
                (-0.5, skycell.pixel_shape[1] - 0.5),
            ]
        ),
        rtol=1e-5,
    )


def test_exclusive_core():
    rng = np.random.default_rng()
    ra = rng.standard_normal(1000) * 360.0
    dec = rng.standard_normal(1000) * 180.0 - 90
    radec = np.stack([ra, dec], axis=1)

    vectorpoints = np.stack(sgv.lonlat_to_vector(radec[:, 0], radec[:, 1]), axis=1)
    if not isinstance(vectorpoints, np.ndarray):
        vectorpoints = np.ndarray([vectorpoints])
    vectorpoints = sgv.normalize_vector(vectorpoints)

    skycells_kdtree = sc.SKYMAP.skycells_kdtree

    point_nearest_skycell_indices = skycells_kdtree.query(
        vectorpoints, k=16, distance_upper_bound=sc.SkyCell.length
    )[1]

    skycells = []
    for point_index, skycell_indices in enumerate(point_nearest_skycell_indices):
        for skycell_index in skycell_indices:
            if skycell_index != len(sc.SKYMAP.model.skycells):
                skycell = sc.SkyCell(skycell_index)
                
                if skycell.core_contains(radec[point_index, :]):
                    skycells.append(skycell)
                    break

    assert len(skycells) == radec.shape[0]
