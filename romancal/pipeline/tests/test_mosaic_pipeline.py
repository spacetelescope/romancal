"""Unit tests for the mosaic pipeline"""
import numpy as np

import romancal.pipeline.mosaic_pipeline as mp

def test_wcsinfo_to_wcs():
    """Test integrity of wcsinfo_to_wcs"""
    wcsinfo = {
        'ra_ref': 269.83219987378925,
        'dec_ref': 66.04081466149024,
        'x_ref': 2069.0914958388985,
        'y_ref': 2194.658767532754,
        'rotation_matrix': [[-0.9999964196507396, -0.00267594575838714], [-0.00267594575838714, 0.9999964196507396]],
        'pixel_scale': 3.036307317109957e-05,
        'pixel_shape': [4389, 4138],
        'ra_center': 269.82284964811464,
        'dec_center': 66.0369888162117,
        'ra_corn1': 269.98694025887136,
        'dec_corn1': 65.97426875366378,
        'ra_corn2': 269.98687579251805,
        'dec_corn2': 66.09988065827382,
        'ra_corn3': 269.6596332616879,
        'dec_corn3': 65.97389321243348,
        'ra_corn4': 269.6579498847431,
        'dec_corn4': 66.099533603104,
        'orientat': 359.8466793994546
    }

    wcs = mp.wcsinfo_to_wcs(wcsinfo)

    assert np.allclose(wcs(wcsinfo['x_ref'], wcsinfo['y_ref']), (wcsinfo['ra_ref'], wcsinfo['dec_ref']))
    assert np.allclose(wcs(4389 / 2., 4138 / 2.), (wcsinfo['ra_center'], wcsinfo['dec_center']))
    assert np.allclose(wcs(0., 0.), (wcsinfo['ra_corn1'], wcsinfo['dec_corn1']))
    assert np.allclose(wcs(0., wcsinfo['pixel_shape'][1] - 1.), (wcsinfo['ra_corn2'], wcsinfo['dec_corn2']))
    assert np.allclose(wcs(wcsinfo['pixel_shape'][0], 0.), (wcsinfo['ra_corn3'], wcsinfo['dec_corn3']))
    assert np.allclose(wcs(wcsinfo['pixel_shape'][0], wcsinfo['pixel_shape'][1]), (wcsinfo['ra_corn4'], wcsinfo['dec_corn4']))
