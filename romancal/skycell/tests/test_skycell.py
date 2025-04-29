"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

import romancal.skycell.skymap as sc


def assert_allclose_lonlat(actual: np.ndarray, desired: np.ndarray, rtol=1e-7, atol=0):
    assert_allclose(
        np.array(actual) % 360, np.array(desired) % 360, rtol=rtol, atol=atol
    )


@pytest.fixture(autouse=True)
def override_skymap(monkeypatch):
    """
    For the tests in this file, monkeypatch the global
    skymap path to a smaller subset to allow these tests
    to run without access to the full skymap from CRDS.
    """
    monkeypatch.setattr(sc.SKYMAP, "path", Path(__file__).parent / "skymap_subset.asdf")
    yield


def test_skycell_init():
    skycell = sc.SkyCell.from_name("225p90x30y51")

    assert skycell == sc.SkyCell(107)

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
        sc.SkyCell.from_name("r274dp63x31y81")

    with pytest.raises(ValueError):
        sc.SkyCell.from_name("notaskycellname")


def test_skycell_from_projregion():
    projregion = sc.ProjectionRegion(0)

    assert projregion.skycells[100] == sc.SkyCell.from_name("225p90x30y44").data

    assert projregion.skycell_indices[-1] != sc.ProjectionRegion(1).skycell_indices[0]


def test_projregion_from_skycell():
    skycell = sc.SkyCell.from_name("225p90x30y51")

    assert skycell.projection_region == sc.ProjectionRegion(0)

    assert sc.ProjectionRegion.from_skycell_index(107) == sc.ProjectionRegion(0)

    assert sc.ProjectionRegion.from_skycell_index(0) == sc.ProjectionRegion(0)


@pytest.mark.parametrize(
    "name",
    [
        "000p86x30y34",
        "000p86x50y65",
        "000p86x59y38",
        "000p86x61y68",
        "045p86x29y34",
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
def test_skycell_wcs_pixel_to_world(name):
    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs

    # forward transform to radec corners
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
def test_skycell_wcs_world_to_pixel(name):
    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs

    # inverse transform to pixel corners
    assert_allclose(
        np.array(wcs.invert(*skycell.radec_corners.T)).T,
        (
            [
                (0.5, 0.5),
                (skycell.pixel_shape[0] + 0.5, 0.5),
                (skycell.pixel_shape[0] + 0.5, skycell.pixel_shape[1] + 0.5),
                (0.5, skycell.pixel_shape[1] + 0.5),
            ]
        ),
        # TODO; figure out why this is so messy in the case of polar caps
        rtol=1e-5 if not skycell.projection_region.is_polar else 1e-1,
    )


@pytest.mark.parametrize(
    "name",
    [
        "000p86x30y34",
        "000p86x50y65",
        "000p86x59y38",
        "000p86x61y68",
        "045p86x29y34",
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
def test_skycell_wcsinfo(name):
    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs
    wcsinfo = skycell.wcsinfo

    assert_allclose_lonlat(
        wcs(wcsinfo["x_ref"], wcsinfo["y_ref"]),
        (wcsinfo["ra_ref"], wcsinfo["dec_ref"]),
        rtol=1e-5,
    )
    assert_allclose_lonlat(
        wcs(
            (wcsinfo["pixel_shape"][0] / 2.0) + 0.5,
            (wcsinfo["pixel_shape"][1] / 2.0) + 0.5,
        ),
        (wcsinfo["ra_center"], wcsinfo["dec_center"]),
        rtol=1e-5,
    )

    assert_allclose_lonlat(
        np.array(
            wcs(
                *np.array(
                    [
                        (0.5, 0.5),
                        (wcsinfo["pixel_shape"][0] + 0.5, 0.5),
                        (
                            wcsinfo["pixel_shape"][0] + 0.5,
                            wcsinfo["pixel_shape"][1] + 0.5,
                        ),
                        (0.5, wcsinfo["pixel_shape"][1] + 0.5),
                    ]
                ).T
            )
        ).T,
        (
            [
                (wcsinfo["ra_corn1"], wcsinfo["dec_corn1"]),
                (wcsinfo["ra_corn2"], wcsinfo["dec_corn2"]),
                (wcsinfo["ra_corn3"], wcsinfo["dec_corn3"]),
                (wcsinfo["ra_corn4"], wcsinfo["dec_corn4"]),
            ]
        ),
        rtol=1e-5,
    )
