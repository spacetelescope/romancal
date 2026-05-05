"""
Tests for the ``available_properties``, ``properties``, and
``requested_properties`` API of the source catalog sub-catalogs.
"""

from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest

from romancal.source_catalog.aperture import ApertureCatalog
from romancal.source_catalog.daofind import DAOFindCatalog
from romancal.source_catalog.neighbors import NNCatalog

NN_AVAILABLE_PROPERTIES = ("nn_label", "nn_distance")
DAO_AVAILABLE_PROPERTIES = ("sharpness", "roundness1")


@pytest.fixture
def nn_inputs():
    label = np.array([1, 2, 3], dtype=np.int32)
    xypos = np.array([[10.0, 10.0], [20.0, 20.0], [30.0, 30.0]])
    return label, xypos, xypos.copy(), 0.1 * u.arcsec


class TestNNCatalogProperties:
    def test_available_properties_is_complete(self):
        assert set(NN_AVAILABLE_PROPERTIES) == set(NNCatalog.available_properties)

    def test_default_properties_match_available(self, nn_inputs):
        cat = NNCatalog(*nn_inputs)
        assert list(cat.properties) == list(NNCatalog.available_properties)

    @pytest.mark.parametrize(
        "requested, expected",
        [
            (["nn_label"], ["nn_label"]),
            (["nn_distance"], ["nn_distance"]),
            (["nn_label", "nn_distance"], ["nn_label", "nn_distance"]),
            (["unrelated"], []),
            ([], []),
        ],
    )
    def test_requested_properties_filter(self, nn_inputs, requested, expected):
        cat = NNCatalog(*nn_inputs, requested_properties=requested)
        assert cat.properties == expected

    def test_properties_preserve_available_order(self, nn_inputs):
        # Request in reversed order — output order must match
        # available_properties
        cat = NNCatalog(*nn_inputs, requested_properties=["nn_distance", "nn_label"])
        assert cat.properties == list(NNCatalog.available_properties)


@pytest.fixture
def daofind_inputs():
    rng = np.random.default_rng(42)
    data = rng.standard_normal((50, 50))
    xypos = np.array([[25.0, 25.0]])
    return data, xypos, 2.0


class TestDAOFindCatalogProperties:
    def test_available_properties_is_complete(self):
        assert set(DAO_AVAILABLE_PROPERTIES) == set(DAOFindCatalog.available_properties)

    def test_default_properties_match_available(self, daofind_inputs):
        cat = DAOFindCatalog(*daofind_inputs)
        assert list(cat.properties) == list(DAOFindCatalog.available_properties)

    @pytest.mark.parametrize(
        "requested, expected",
        [
            (["sharpness"], ["sharpness"]),
            (["roundness1"], ["roundness1"]),
            (["sharpness", "roundness1"], ["sharpness", "roundness1"]),
            (["nn_label"], []),
            ([], []),
        ],
    )
    def test_requested_properties_filter(self, daofind_inputs, requested, expected):
        cat = DAOFindCatalog(*daofind_inputs, requested_properties=requested)
        assert cat.properties == expected


def _make_aperture_inputs():
    """
    Build a minimal model + xypos for ApertureCatalog.
    """
    data = np.zeros((50, 50)) << u.nJy
    err = np.ones((50, 50)) << u.nJy
    model = SimpleNamespace(data=data, err=err)
    pixel_scale = 0.11 * u.arcsec
    xypos = np.array([[25.0, 25.0]])
    return model, xypos, pixel_scale


class TestApertureCatalogProperties:
    def test_static_aperture_flux_colnames(self):
        # The flux column names are deterministic from the radii in arcsec
        names = ApertureCatalog.aperture_flux_colnames_for_radii()
        assert names == [
            "aper01_flux",
            "aper02_flux",
            "aper04_flux",
            "aper08_flux",
            "aper16_flux",
        ]

    def test_available_properties_includes_flux_and_err_and_bkg(self):
        model, xypos, pixel_scale = _make_aperture_inputs()
        cat = ApertureCatalog(model, xypos, pixel_scale)
        available = set(cat.available_properties)
        # Every flux column should appear with a matching `_err`
        for name in ApertureCatalog.aperture_flux_colnames_for_radii():
            assert name in available
            assert f"{name}_err" in available
        assert "aper_bkg_flux" in available
        assert "aper_bkg_flux_err" in available

    def test_default_properties_match_available(self):
        model, xypos, pixel_scale = _make_aperture_inputs()
        cat = ApertureCatalog(model, xypos, pixel_scale)
        assert list(cat.properties) == list(cat.available_properties)

    def test_requested_properties_filter(self):
        model, xypos, pixel_scale = _make_aperture_inputs()
        requested = ["aper02_flux", "aper02_flux_err", "aper_bkg_flux"]
        cat = ApertureCatalog(model, xypos, pixel_scale, requested_properties=requested)
        assert set(cat.properties) == set(requested)
        # Order is preserved from `available_properties`
        assert cat.properties == [n for n in cat.available_properties if n in requested]

    def test_requested_properties_unknown_ignored(self):
        model, xypos, pixel_scale = _make_aperture_inputs()
        cat = ApertureCatalog(
            model, xypos, pixel_scale, requested_properties=["unrelated"]
        )
        assert cat.properties == []

    def test_is_extended_no_ee_spline_returns_all_false(self):
        model, xypos, pixel_scale = _make_aperture_inputs()
        cat = ApertureCatalog(model, xypos, pixel_scale)
        result = cat.is_extended
        assert result.dtype == bool
        assert result.shape == (xypos.shape[0],)
        assert not result.any()
