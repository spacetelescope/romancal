"""Unit tests for skycell wcs functions"""

import os
from pathlib import Path

import numpy as np
import pytest

os.environ["SKYMAP_PATH"] = str(Path(__file__).parent / "skymap_subset.asdf")
from romancal.skycell.skymap import SkyCell, wcsinfo_to_gwcs


def test_skycell_index():
    skycell = SkyCell.from_data(
        np.void(
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
    )

    assert skycell.data == SkyCell(107).data
    assert skycell.data == SkyCell.from_name("225p90x30y51").data

    with pytest.raises(ValueError):
        SkyCell.from_name("r274dp63x31y81")

    with pytest.raises(ValueError):
        SkyCell.from_name("notaskycellname")


def test_skycell_to_wcs():
    """Test integrity of skycell.wcs"""

    skycell = SkyCell.from_name("225p90x30y51")

    assert np.allclose(
        skycell.wcs.wcs_pix2world([(0.0, 0.0)], 0),
        skycell.radec_corners[0],
    )
    assert np.allclose(
        skycell.wcs.wcs_pix2world([(skycell.pixel_shape[0] - 1, 0.0)], 0),
        skycell.radec_corners[1],
    )
    assert np.allclose(
        skycell.wcs.wcs_pix2world(
            [(skycell.pixel_shape[0] - 1, skycell.pixel_shape[1] - 1)], 0
        ),
        skycell.radec_corners[2],
    )
    assert np.allclose(
        skycell.wcs.wcs_pix2world([(0.0, skycell.pixel_shape[1] - 1)], 0),
        skycell.radec_corners[3],
    )


def test_wcsinfo_to_wcs():
    """Test integrity of wcsinfo_to_wcs"""
    wcsinfo = {
        "ra_ref": 269.83219987378925,
        "dec_ref": 66.04081466149024,
        "x_ref": 2069.0914958388985,
        "y_ref": 2194.658767532754,
        "rotation_matrix": [
            [-0.9999964196507396, -0.00267594575838714],
            [-0.00267594575838714, 0.9999964196507396],
        ],
        "pixel_scale": 3.036307317109957e-05,
        "pixel_shape": [4389, 4138],
        "ra_center": 269.82284964811464,
        "dec_center": 66.0369888162117,
        "ra_corn1": 269.98694025887136,
        "dec_corn1": 65.97426875366378,
        "ra_corn2": 269.98687579251805,
        "dec_corn2": 66.09988065827382,
        "ra_corn3": 269.6579498847431,
        "dec_corn3": 66.099533603104,
        "ra_corn4": 269.6596332616879,
        "dec_corn4": 65.97389321243348,
        "orientat": 359.8466793994546,
    }

    wcs = wcsinfo_to_gwcs(wcsinfo)

    assert np.allclose(
        wcs(wcsinfo["x_ref"], wcsinfo["y_ref"]), (wcsinfo["ra_ref"], wcsinfo["dec_ref"])
    )
    assert np.allclose(
        wcs(4389 / 2.0, 4138 / 2.0), (wcsinfo["ra_center"], wcsinfo["dec_center"])
    )
    assert np.allclose(wcs(0.0, 0.0), (wcsinfo["ra_corn1"], wcsinfo["dec_corn1"]))
    assert np.allclose(
        wcs(0.0, wcsinfo["pixel_shape"][1]), (wcsinfo["ra_corn2"], wcsinfo["dec_corn2"])
    )
    assert np.allclose(
        wcs(wcsinfo["pixel_shape"][0], wcsinfo["pixel_shape"][1]),
        (wcsinfo["ra_corn3"], wcsinfo["dec_corn3"]),
    )
    assert np.allclose(
        wcs(wcsinfo["pixel_shape"][0], 0.0), (wcsinfo["ra_corn4"], wcsinfo["dec_corn4"])
    )
