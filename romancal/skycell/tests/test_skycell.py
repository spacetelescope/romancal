"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from romancal.skycell import skymap

DATA_DIRECTORY = Path(__file__).parent / "data"


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(
        np.array(actual) % 360, np.array(desired) % 360, rtol=rtol, atol=atol
    )


@pytest.fixture()
def skymap_subset() -> skymap.SkyMap:
    """
    smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    return skymap.SkyMap(DATA_DIRECTORY / "skymap_subset.asdf")


def test_skycell_from_name(skymap_subset):
    skycell = skymap.SkyCell.from_name("135p90x50y57", skymap=skymap_subset)

    assert skycell == skymap.SkyCell(999, skymap=skymap_subset)

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

    with pytest.raises(ValueError):
        skymap.SkyCell.from_name("r274dp63x63y81", skymap=skymap_subset)

    with pytest.raises(ValueError):
        skymap.SkyCell.from_name("notaskycellname", skymap=skymap_subset)


def test_skycell_from_asn(skymap_subset):
    skycell = skymap.SkyCell.from_asn(
        DATA_DIRECTORY / "L3_mosaic_asn.json", skymap=skymap_subset
    )
    assert skycell.name == "000p86x69y62"

    with pytest.raises(ValueError):
        skymap.SkyCell.from_asn(
            DATA_DIRECTORY / "L3_regtest_asn.json", skymap=skymap_subset
        )

    with pytest.raises(ValueError):
        skymap.SkyCell.from_asn(
            DATA_DIRECTORY / "L3_skycell_mbcat_asn.json", skymap=skymap_subset
        )


def test_skycell_from_projregion(skymap_subset):
    projregion = skymap.ProjectionRegion(0, skymap=skymap_subset)

    assert (
        projregion.skycells[100]
        == skymap.SkyCell.from_name("135p90x30y44", skymap=skymap_subset).data
    )

    assert (
        projregion.skycell_indices[-1]
        != skymap.ProjectionRegion(1, skymap=skymap_subset).skycell_indices[0]
    )


def test_projregion_from_skycell(skymap_subset):
    skycell = skymap.SkyCell.from_name("135p90x50y57", skymap=skymap_subset)

    projregion0 = skymap.ProjectionRegion(0, skymap=skymap_subset)
    projregion1 = skymap.ProjectionRegion(1, skymap=skymap_subset)

    assert skycell.projection_region == projregion0 # this calls CRDS!

    assert (
        skymap.ProjectionRegion.from_skycell_index(107, skymap=skymap_subset) == projregion0
    )

    assert (
        skymap.ProjectionRegion.from_skycell_index(0, skymap=skymap_subset) == projregion0
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


@pytest.mark.parametrize(
    "name",
    [
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
    ],
)
def test_skycell_wcs_pixel_to_world(name, skymap_subset):
    skycell = skymap.SkyCell.from_name(name, skymap=skymap_subset)

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
    ],
)
def test_skycell_wcs_world_to_pixel(name, skymap_subset):
    skycell = skymap.SkyCell.from_name(name, skymap=skymap_subset)

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


@pytest.mark.parametrize(
    "name",
    [
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
    ],
)
def test_skycell_wcsinfo(name, skymap_subset):
    skycell = skymap.SkyCell.from_name(name, skymap=skymap_subset)

    wcs = skycell.wcs
    wcsinfo = skycell.wcs_info

    assert_allclose_lonlat(
        wcs(
            (wcsinfo["nx"] / 2.0) - 0.5,
            (wcsinfo["ny"] / 2.0) - 0.5,
        ),
        (wcsinfo["ra_center"], wcsinfo["dec_center"]),
        rtol=1e-7,
    )

    assert_allclose_lonlat(
        np.array(
            wcs(
                *np.array(
                    [
                        (-0.5, -0.5),
                        (wcsinfo["nx"] - 0.5, -0.5),
                        (
                            wcsinfo["nx"] - 0.5,
                            wcsinfo["ny"] - 0.5,
                        ),
                        (-0.5, wcsinfo["ny"] - 0.5),
                    ]
                ).T,
                with_bounding_box=False,
            )
        ).T,
        skycell.radec_corners,
        rtol=1e-7,
    )
