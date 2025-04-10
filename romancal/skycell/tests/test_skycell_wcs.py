"""Unit tests for skycell wcs functions"""

import numpy as np

from romancal.skycell.match import skycell_to_wcs, wcsinfo_to_wcs


def test_skycell_to_wcs():
    """Test integrity of skycell_to_wcs"""

    skycell = np.void(
        (
            "r274dp63x31y81",
            269.7783307416819,
            66.04965143695566,
            1781.5,
            1781.5,
            355.9788,
            3564,
            3564,
            67715.5,
            -110484.5,
            269.6657957545588,
            65.9968687812357,
            269.6483032937494,
            66.09523979539262,
            269.89132874168854,
            66.10234971630734,
            269.9079118635897,
            66.00394719483091,
            0.1,
            274.2857142857143,
            63.0,
            0.0,
            463181,
        ),
        dtype=[
            ("name", "<U20"),
            ("ra_center", "<f8"),
            ("dec_center", "<f8"),
            ("x_center", "<f4"),
            ("y_center", "<f4"),
            ("orientat", "<f4"),
            ("nx", "<i4"),
            ("ny", "<i4"),
            ("x0_projection", "<f4"),
            ("y0_projection", "<f4"),
            ("ra_corn1", "<f8"),
            ("dec_corn1", "<f8"),
            ("ra_corn2", "<f8"),
            ("dec_corn2", "<f8"),
            ("ra_corn3", "<f8"),
            ("dec_corn3", "<f8"),
            ("ra_corn4", "<f8"),
            ("dec_corn4", "<f8"),
            ("pixel_scale", "<f4"),
            ("ra_projection_center", "<f8"),
            ("dec_projection_center", "<f8"),
            ("orientat_projection_center", "<f4"),
            ("index", "<i8"),
        ],
    )

    wcs = skycell_to_wcs(skycell)

    assert np.allclose(
        wcs(
            skycell["x0_projection"], skycell["y0_projection"], with_bounding_box=False
        ),
        (skycell["ra_projection_center"], skycell["dec_projection_center"]),
    )
    assert np.allclose(
        wcs(skycell["x_center"], skycell["y_center"]),
        (skycell["ra_center"], skycell["dec_center"]),
    )
    assert np.allclose(wcs(0.0, 0.0), (skycell["ra_corn1"], skycell["dec_corn1"]))
    assert np.allclose(
        wcs(0.0, skycell["ny"] - 1), (skycell["ra_corn2"], skycell["dec_corn2"])
    )
    assert np.allclose(
        wcs(skycell["nx"] - 1, skycell["ny"] - 1),
        (skycell["ra_corn3"], skycell["dec_corn3"]),
    )
    assert np.allclose(
        wcs(skycell["nx"] - 1, 0.0), (skycell["ra_corn4"], skycell["dec_corn4"])
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
