"""Unit tests for skycell wcs functions"""

import numpy as np
import pytest

from romancal.skycell.skycells import SkyCell, wcsinfo_to_wcs


def test_skycell_index():
    skycell = SkyCell.from_data(
        np.void(
            (
                "206p65x36y36",
                203.4088863995,
                63.5573382312,
                2.3053994,
                69699.5,
                69699.5,
                203.3264435081,
                63.5177894396,
                203.4975627644,
                63.5205352377,
                203.4915622721,
                63.5968414781,
                203.3199808485,
                63.5940863621,
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

    assert skycell.data == SkyCell(400000).data
    assert skycell.data == SkyCell.from_name("206p65x36y36").data

    with pytest.raises(ValueError):
        SkyCell.from_name("r274dp63x31y81")

    with pytest.raises(ValueError):
        SkyCell.from_name("notaskycellname")


def test_skycell_to_wcs():
    """Test integrity of skycell.wcs"""

    # skycell = SkyCell.from_data(
    #     np.void(
    #         (
    #             "206p65x36y36",
    #             203.4088863995,
    #             63.5573382312,
    #             2.3053994,
    #             69699.5,
    #             69699.5,
    #             203.3264435081,
    #             63.5177894396,
    #             203.4975627644,
    #             63.5205352377,
    #             203.4915622721,
    #             63.5968414781,
    #             203.3199808485,
    #             63.5940863621,
    #         ),
    #         dtype=[
    #             ("name", "<U16"),
    #             ("ra_center", "<f8"),
    #             ("dec_center", "<f8"),
    #             ("orientat", "<f4"),
    #             ("x_tangent", "<f8"),
    #             ("y_tangent", "<f8"),
    #             ("ra_corn1", "<f8"),
    #             ("dec_corn1", "<f8"),
    #             ("ra_corn2", "<f8"),
    #             ("dec_corn2", "<f8"),
    #             ("ra_corn3", "<f8"),
    #             ("dec_corn3", "<f8"),
    #             ("ra_corn4", "<f8"),
    #             ("dec_corn4", "<f8"),
    #         ],
    #     )
    # )

    skycell = SkyCell.from_name("206p65x36y36")

    assert np.allclose(
        skycell.wcs(0.0, skycell.pixel_shape[1] - 1),
        skycell.radec_corners[0],
    )
    assert np.allclose(
        skycell.wcs(skycell.pixel_shape[0] - 1, skycell.pixel_shape[1] - 1),
        skycell.radec_corners[1],
    )
    assert np.allclose(
        skycell.wcs(skycell.pixel_shape[0] - 1, skycell.pixel_shape[1] - 1),
        skycell.radec_corners[2],
    )
    assert np.allclose(skycell.wcs(0.0, 0.0), skycell.radec_corners[3])


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

    wcs = wcsinfo_to_wcs(wcsinfo)

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
