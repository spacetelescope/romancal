"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

import romancal.skycell.skymap as sc

DATA_DIRECTORY = Path(__file__).parent / "data"


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(
        np.array(actual) % 360, np.array(desired) % 360, rtol=rtol, atol=atol
    )


@pytest.fixture()
def skymap_subset() -> sc.SkyMap:
    """
    smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    return sc.SkyMap(DATA_DIRECTORY / "skymap_subset.asdf")


def test_skycell_from_name(skymap_subset):
    skycell = sc.SkyCell.from_name("225p90x30y51", skymap=skymap_subset)

    assert skycell == sc.SkyCell(107, skymap=skymap_subset)

    assert skycell.data == np.void(
        (
            "225p90x30y51",
            312.1369696579,
            88.5318370582,
            -87.13697,
            98499.5,
            -2300.5,
            313.6617771309,
            88.4950897773,
            313.5902417923,
            88.5714066707,
            310.5350393365,
            88.5674932367,
            310.7608411964,
            88.4913746368,
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
        sc.SkyCell.from_name("r274dp63x63y81", skymap=skymap_subset)

    with pytest.raises(ValueError):
        sc.SkyCell.from_name("notaskycellname", skymap=skymap_subset)


def test_skycell_from_asn(skymap_subset):
    skycell = sc.SkyCell.from_asn(
        DATA_DIRECTORY / "L3_mosaic_asn.json", skymap=skymap_subset
    )
    assert skycell.name == "225p90x49y67"

    with pytest.raises(ValueError):
        sc.SkyCell.from_asn(
            DATA_DIRECTORY / "L3_regtest_asn.json", skymap=skymap_subset
        )

    with pytest.raises(ValueError):
        sc.SkyCell.from_asn(
            DATA_DIRECTORY / "L3_skycell_mbcat_asn.json", skymap=skymap_subset
        )


def test_skycell_from_projregion(skymap_subset):
    projregion = sc.ProjectionRegion(0, skymap=skymap_subset)

    assert (
        projregion.skycells[100]
        == sc.SkyCell.from_name("225p90x30y44", skymap=skymap_subset).data
    )

    assert (
        projregion.skycell_indices[-1]
        != sc.ProjectionRegion(1, skymap=skymap_subset).skycell_indices[0]
    )


def test_projregion_from_skycell(skymap_subset):
    skycell = sc.SkyCell.from_name("225p90x30y51")

    projregion0 = sc.ProjectionRegion(0, skymap=skymap_subset)
    projregion1 = sc.ProjectionRegion(1, skymap=skymap_subset)

    assert skycell.projection_region == projregion0

    assert (
        sc.ProjectionRegion.from_skycell_index(107, skymap=skymap_subset) == projregion0
    )

    assert (
        sc.ProjectionRegion.from_skycell_index(0, skymap=skymap_subset) == projregion0
    )

    assert (
        sc.ProjectionRegion.from_skycell_index(
            projregion0.data["skycell_end"], skymap=skymap_subset
        )
        == projregion1
    )

    with pytest.raises(KeyError):
        sc.ProjectionRegion.from_skycell_index(-1, skymap=skymap_subset)

    with pytest.raises(KeyError):
        sc.ProjectionRegion.from_skycell_index(10000, skymap=skymap_subset)


@pytest.mark.parametrize(
    "name",
    [
        "000p86x30y34",
        "000p86x50y65",
        "000p86x59y38",
        "000p86x61y68",
        "045p86x29y34",
        # north pole
        "225p90x25y49",
        "225p90x30y51",
        "225p90x33y62",
        "225p90x39y33",
        "225p90x43y65",
        "225p90x48y41",
        "225p90x52y59",
        "225p90x57y35",
        "225p90x61y67",
        "225p90x67y38",
    ],
)
def test_skycell_wcs_pixel_to_world(name, skymap_subset):
    skycell = sc.SkyCell.from_name(name, skymap=skymap_subset)

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
        "000p86x61y68",
        "045p86x29y34",
        # north pole
        "225p90x25y49",
        "225p90x30y51",
        "225p90x33y62",
        "225p90x39y33",
        "225p90x43y65",
        "225p90x48y41",
        "225p90x52y59",
        "225p90x57y35",
        "225p90x61y67",
        "225p90x67y38",
    ],
)
def test_skycell_wcs_world_to_pixel(name, skymap_subset):
    skycell = sc.SkyCell.from_name(name, skymap=skymap_subset)

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


@pytest.mark.parametrize(
    "name",
    [
        "000p86x30y34",
        "000p86x50y65",
        "000p86x59y38",
        "000p86x61y68",
        "045p86x29y34",
        # north pole
        "225p90x25y49",
        "225p90x30y51",
        "225p90x33y62",
        "225p90x39y33",
        "225p90x43y65",
        "225p90x48y41",
        "225p90x52y59",
        "225p90x57y35",
        "225p90x61y67",
        "225p90x67y38",
    ],
)
def test_skycell_wcsinfo(name, skymap_subset):
    skycell = sc.SkyCell.from_name(name, skymap=skymap_subset)

    wcs = skycell.wcs
    wcsinfo = skycell.wcs_info

    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose_lonlat(
        wcs(
            (wcsinfo["nx"] / 2.0) + 0.5,
            (wcsinfo["ny"] / 2.0) + 0.5,
        ),
        (wcsinfo["ra_center"], wcsinfo["dec_center"]),
        rtol=1e-7,
    )

    # TODO: the corners in the reference file currently use FITS convention (pixel + 0.5) instead of (pixel - 0.5)
    assert_allclose_lonlat(
        np.array(
            wcs(
                *np.array(
                    [
                        (0.5, 0.5),
                        (wcsinfo["nx"] + 0.5, 0.5),
                        (
                            wcsinfo["nx"] + 0.5,
                            wcsinfo["ny"] + 0.5,
                        ),
                        (0.5, wcsinfo["ny"] + 0.5),
                    ]
                ).T,
                with_bounding_box=False,
            )
        ).T,
        skycell.radec_corners,
        rtol=1e-7,
    )
