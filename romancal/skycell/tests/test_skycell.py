"""Unit tests for skycell functions"""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

import romancal.skycell.skymap as sc


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
    ["000p86x30y34", "000p86x50y65", "000p86x61y68", "045p86x29y34", "225p90x30y51"],
)
def test_skycell_to_wcs(name):
    """Test integrity of skycell.wcs"""

    skycell = sc.SkyCell.from_name(name)

    wcs = skycell.wcs

    # forward transform to radec corners
    assert_allclose(
        np.array(
            [
                wcs(0.0, 0.0),
                wcs(skycell.pixel_shape[0] - 1, 0.0),
                wcs(skycell.pixel_shape[0] - 1, skycell.pixel_shape[1] - 1),
                wcs(0.0, skycell.pixel_shape[1] - 1),
            ]
        ),
        skycell.radec_corners,
    )

    # inverse transform to pixel corners
    assert_allclose(
        np.array([wcs.invert(*corner) for corner in skycell.radec_corners]),
        np.array(
            [
                (0.0, 0.0),
                (skycell.pixel_shape[0] - 1, 0.0),
                (skycell.pixel_shape[0] - 1, skycell.pixel_shape[1] - 1),
                (0.0, skycell.pixel_shape[1] - 1),
            ]
        ),
    )

    wcsinfo = skycell.wcsinfo

    assert_allclose(
        wcs(wcsinfo["x_ref"], wcsinfo["y_ref"]), (wcsinfo["ra_ref"], wcsinfo["dec_ref"])
    )
    assert_allclose(
        wcs(wcsinfo["pixel_shape"][0] / 2.0, wcsinfo["pixel_shape"][1] / 2.0),
        (wcsinfo["ra_center"], wcsinfo["dec_center"]),
    )
    assert_allclose(wcs(0.0, 0.0), (wcsinfo["ra_corn1"], wcsinfo["dec_corn1"]))
    assert_allclose(
        wcs(0.0, wcsinfo["pixel_shape"][1]), (wcsinfo["ra_corn2"], wcsinfo["dec_corn2"])
    )
    assert_allclose(
        wcs(wcsinfo["pixel_shape"][0], wcsinfo["pixel_shape"][1]),
        (wcsinfo["ra_corn3"], wcsinfo["dec_corn3"]),
    )
    assert_allclose(
        wcs(wcsinfo["pixel_shape"][0], 0.0), (wcsinfo["ra_corn4"], wcsinfo["dec_corn4"])
    )
