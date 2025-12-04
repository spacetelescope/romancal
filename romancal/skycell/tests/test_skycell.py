"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from romancal.skycell import skymap

DATA_DIRECTORY = Path(__file__).parent / "data"

TEST_SKYCELLS = [
    "000p86x30y34",
    "000p86x50y65",
    "000p86x59y38",
    # north pole
    "135p90x25y49",
    "135p90x30y51",
    "135p90x33y62",
    "135p90x39y33",
    "135p90x43y65",
    "135p90x48y41",
    "135p90x52y59",
    "135p90x57y35",
    "135p90x61y67",
    "135p90x67y38",
]


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(actual, desired)


@pytest.fixture()
def skymap_subset() -> skymap.SkyMap:
    """
    smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    return skymap.SkyMap(DATA_DIRECTORY / "skymap_subset.asdf")


def test_skycell_deprecation(skymap_subset):
    with pytest.warns(expected_warning=DeprecationWarning):
        skymap.SkyCell.from_name("135p90x50y57", skymap=skymap_subset)


def test_skycells(skymap_subset):
    skycells = skymap.SkyCells.from_names(TEST_SKYCELLS, skymap=skymap_subset)

    assert sorted(skycells.names) == sorted(TEST_SKYCELLS)

    assert skycells.radec_corners.shape == (len(TEST_SKYCELLS), 4, 2)
    assert skycells.vectorpoint_corners.shape == (len(TEST_SKYCELLS), 4, 3)

    assert skycells.radec_centers.shape == (len(TEST_SKYCELLS), 2)
    assert skycells.vectorpoint_centers.shape == (len(TEST_SKYCELLS), 3)

    assert len(skycells.polygons) == len(TEST_SKYCELLS)


def test_skycell_from_name(skymap_subset):
    skycell = skymap.SkyCells.from_names(["135p90x50y57"], skymap=skymap_subset)

    assert skycell == skymap.SkyCells([999], skymap=skymap_subset)

    assert skycell.data == np.void(
        (
            "135p90x50y57",
            224.99999999999997,
            89.48668040110688,
            -90.0,
            -31100.5,
            2499.5,
            220.4041123081352,
            89.52333943360831,
            221.0384731174898,
            89.44716843824965,
            228.96152688251019,
            89.44716843824965,
            229.59588769186476,
            89.52333943360831,
        ),
        dtype=[
            ("name", "<U16"),
            ("ra_center", "<f8"),
            ("dec_center", "<f8"),
            ("orientat", "<f4"),
            ("x_tangent", "<f8"),
            ("y_tangent", "<f8"),
            ("ra_corn1", "<f8"),
            ("dec_corn1", "<f8"),
            ("ra_corn2", "<f8"),
            ("dec_corn2", "<f8"),
            ("ra_corn3", "<f8"),
            ("dec_corn3", "<f8"),
            ("ra_corn4", "<f8"),
            ("dec_corn4", "<f8"),
        ],
    )

    with pytest.raises(KeyError):
        skymap.SkyCells.from_names(["r274dp63x63y81"], skymap=skymap_subset)

    with pytest.raises(KeyError):
        skymap.SkyCells.from_names(["notaskycellname"], skymap=skymap_subset)


def test_skycell_from_asn(skymap_subset):
    skycell = skymap.SkyCells.from_asns(
        [DATA_DIRECTORY / "L3_mosaic_asn.json"], skymap=skymap_subset
    )
    assert skycell.names == ["000p86x69y62"]

    with pytest.raises(ValueError):
        skymap.SkyCells.from_asns(
            [DATA_DIRECTORY / "L3_regtest_asn.json"], skymap=skymap_subset
        )

    with pytest.raises(ValueError):
        skymap.SkyCells.from_asns(
            [DATA_DIRECTORY / "L3_skycell_mbcat_asn.json"], skymap=skymap_subset
        )
    with pytest.raises(ValueError):
        skymap.SkyCells.from_asns(
            DATA_DIRECTORY.glob("*_asn.json"), skymap=skymap_subset
        )


def test_skycell_from_projregion(skymap_subset):
    projregion = skymap.ProjectionRegion(0, skymap=skymap_subset)

    assert (
        projregion.skycells[100]
        == skymap.SkyCells.from_names(["135p90x30y44"], skymap=skymap_subset).data
    )

    assert (
        projregion.skycell_indices[-1]
        != skymap.ProjectionRegion(1, skymap=skymap_subset).skycell_indices[0]
    )


def test_projregion_from_skycell(skymap_subset):
    skycell = skymap.SkyCells.from_names(["135p90x50y57"], skymap=skymap_subset)

    projregion0 = skymap.ProjectionRegion(0, skymap=skymap_subset)
    projregion1 = skymap.ProjectionRegion(1, skymap=skymap_subset)

    assert len(skycell.projection_regions) == 1
    assert skycell.projection_regions[0] == projregion0.index  # this calls CRDS!

    assert (
        skymap.ProjectionRegion.from_skycell_index(107, skymap=skymap_subset)
        == projregion0
    )

    assert (
        skymap.ProjectionRegion.from_skycell_index(0, skymap=skymap_subset)
        == projregion0
    )

    assert (
        skymap.ProjectionRegion.from_skycell_index(
            projregion0.data["skycell_end"], skymap=skymap_subset
        )
        == projregion1
    )

    with pytest.raises(KeyError):
        skymap.ProjectionRegion.from_skycell_index(-1, skymap=skymap_subset)

    with pytest.raises(KeyError):
        skymap.ProjectionRegion.from_skycell_index(10000, skymap=skymap_subset)


@pytest.mark.parametrize("name", TEST_SKYCELLS)
def test_skycell_wcs_pixel_to_world(name, skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=skymap_subset)

    wcsobj = skycell.wcs[0]

    # forward transform to radec corners
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
def test_skycell_wcs_world_to_pixel(name, skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=skymap_subset)

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


@pytest.mark.parametrize("name", TEST_SKYCELLS)
def test_skycell_wcsinfo(name, skymap_subset):
    skycell = skymap.SkyCells.from_names([name], skymap=skymap_subset)

    wcsobj = skycell.wcs[0]
    wcs_info = skycell.wcs_infos[0]

    assert_allclose(
        wcsobj(
            (wcs_info["nx"] / 2.0) - 0.5,
            (wcs_info["ny"] / 2.0) - 0.5,
        ),
        (wcs_info["ra_center"], wcs_info["dec_center"]),
        rtol=1e-7,
    )

    assert_allclose_lonlat(
        np.array(
            wcsobj(
                *np.array(
                    [
                        (-0.5, -0.5),
                        (wcs_info["nx"] - 0.5, -0.5),
                        (
                            wcs_info["nx"] - 0.5,
                            wcs_info["ny"] - 0.5,
                        ),
                        (-0.5, wcs_info["ny"] - 0.5),
                    ]
                ).T,
                with_bounding_box=False,
            )
        ).T,
        skycell.radec_corners[0],
        rtol=1e-7,
    )
