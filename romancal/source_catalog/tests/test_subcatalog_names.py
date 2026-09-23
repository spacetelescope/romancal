"""
Tests for the ``available_properties``, ``properties``, and
``requested_properties`` API of the source catalog sub-catalogs.
"""

from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.segmentation import SegmentationImage, SourceCatalog

from romancal.source_catalog._aperture import ApertureCatalog
from romancal.source_catalog._daofind import DAOFindCatalog
from romancal.source_catalog._neighbors import NNCatalog

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


def _make_segment_img(shape, xypos):
    """
    Build a segmentation image with a single-pixel segment at each
    source position, labeled 1, 2, ... in order.
    """
    segm = np.zeros(shape, dtype=np.int32)
    for label, (x, y) in enumerate(xypos.astype(int), start=1):
        segm[y, x] = label
    return SegmentationImage(segm), np.arange(1, len(xypos) + 1)


def _make_aperture_inputs():
    """
    Build a minimal model, xypos, pixel area map, segmentation image,
    and labels for ApertureCatalog.
    """
    data = np.zeros((50, 50)) << u.nJy
    err = np.ones((50, 50)) << u.nJy
    model = SimpleNamespace(data=data, err=err)
    pixel_scale = 0.11 * u.arcsec
    xypos = np.array([[25.0, 25.0]])
    pixel_area_map = np.full((50, 50), (pixel_scale**2).to_value(u.sr)) << u.sr
    segm, labels = _make_segment_img(data.shape, xypos)
    return model, xypos, pixel_area_map, segm, labels


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
        model, xypos, area_map, segm, labels = _make_aperture_inputs()
        cat = ApertureCatalog(model, xypos, area_map, segm, labels)
        available = set(cat.available_properties)
        # Every flux column should appear with a matching `_err`
        for name in ApertureCatalog.aperture_flux_colnames_for_radii():
            assert name in available
            assert f"{name}_err" in available
        assert "aper_bkg_flux" in available
        assert "aper_bkg_flux_err" in available

    def test_default_properties_match_available(self):
        model, xypos, area_map, segm, labels = _make_aperture_inputs()
        cat = ApertureCatalog(model, xypos, area_map, segm, labels)
        assert list(cat.properties) == list(cat.available_properties)

    def test_requested_properties_filter(self):
        model, xypos, area_map, segm, labels = _make_aperture_inputs()
        requested = ["aper02_flux", "aper02_flux_err", "aper_bkg_flux"]
        cat = ApertureCatalog(
            model, xypos, area_map, segm, labels, requested_properties=requested
        )
        assert set(cat.properties) == set(requested)
        # Order is preserved from `available_properties`
        assert cat.properties == [n for n in cat.available_properties if n in requested]

    def test_requested_properties_unknown_ignored(self):
        model, xypos, area_map, segm, labels = _make_aperture_inputs()
        cat = ApertureCatalog(
            model, xypos, area_map, segm, labels, requested_properties=["unrelated"]
        )
        assert cat.properties == []

    def test_is_extended_no_ee_spline_returns_all_false(self):
        model, xypos, area_map, segm, labels = _make_aperture_inputs()
        cat = ApertureCatalog(model, xypos, area_map, segm, labels)
        result = cat.is_extended
        assert result.dtype == bool
        assert result.shape == (xypos.shape[0],)
        assert not result.any()


class TestApertureRadiiVaryWithPixelArea:
    """
    The aperture radii must track the local pixel solid angle so that
    every source is measured through the same sky aperture.
    """

    def _catalog_with_area_gradient(self):
        data = np.zeros((50, 50)) << u.nJy
        err = np.ones((50, 50)) << u.nJy
        model = SimpleNamespace(data=data, err=err)
        pixel_scale = 0.11 * u.arcsec
        nominal = (pixel_scale**2).to_value(u.sr)

        # Pixel area increases by 10% from left to right
        area_map = np.empty((50, 50))
        area_map[:] = nominal * np.linspace(0.95, 1.05, 50)[np.newaxis, :]
        area_map = area_map << u.sr

        xypos = np.array([[5.0, 25.0], [25.0, 25.0], [45.0, 25.0]])
        segm, labels = _make_segment_img(data.shape, xypos)
        return ApertureCatalog(model, xypos, area_map, segm, labels)

    def test_radii_scale_as_inverse_sqrt_area(self):
        cat = self._catalog_with_area_gradient()
        radii_pix = cat.aperture_radii["circle_pix"]
        areas = cat._source_pixel_area

        assert radii_pix.shape == (len(cat.aperture_radii["circle"]), 3)

        # A source on a larger pixel needs fewer pixels of radius
        assert np.all(np.diff(radii_pix, axis=1) < 0)

        # The sky radius implied by each pixel radius must be constant
        sky_radii = radii_pix * np.sqrt(areas.to_value(u.arcsec**2))[np.newaxis, :]
        expected = cat.aperture_radii["circle"].to_value(u.arcsec)[:, np.newaxis]
        assert np.allclose(sky_radii, expected, rtol=1e-6)


def test_subtract_local_bkg_recovers_source_flux():
    """
    Subtracting the local background must remove a uniform background
    and leave only the source flux, even where the pixel area varies.
    """
    shape = (60, 120)
    bkg_per_pixel = 10.0
    source_flux = 50.0

    # Pixel area varies by 4% across the image, so the sources fall in
    # different radius bins
    nominal = ((0.11 * u.arcsec) ** 2).to_value(u.sr)
    area_map = nominal * np.linspace(0.98, 1.02, shape[1])[np.newaxis, :]
    area_map = np.broadcast_to(area_map, shape) << u.sr

    # Separated by more than the outer annulus radius (~25 pixels)
    xypos = np.array([[30.0, 30.0], [60.0, 30.0], [90.0, 30.0]])
    data = np.full(shape, bkg_per_pixel)
    for x, y in xypos.astype(int):
        data[y, x] += source_flux
    model = SimpleNamespace(data=data << u.nJy, err=np.ones(shape) << u.nJy)
    segm, labels = _make_segment_img(shape, xypos)

    cat = ApertureCatalog(model, xypos, area_map, segm, labels)
    assert len(cat._radius_bins) > 1
    cat.calc_aperture_photometry(subtract_local_bkg=True)

    for name in cat.aperture_flux_colnames:
        assert u.allclose(getattr(cat, name), source_flux * u.nJy, atol=1e-2 * u.nJy)


class TestApertureNeighborMasking:
    """
    Pixels belonging to neighboring sources in the segmentation
    image must be excluded from the circular apertures, matching the
    ``aperture_mask_method="mask"`` used by the segmentation catalog.
    Pixels belonging to any source, including the target, must be
    excluded from the background annulus.
    """

    pixel_scale = 0.11 * u.arcsec
    shape = (80, 80)

    def _uniform_area_map(self):
        area = (self.pixel_scale**2).to_value(u.sr)
        return np.full(self.shape, area) << u.sr

    def test_neighbor_pixels_excluded_from_aperture_flux(self):
        # The target is a single zero-valued pixel. A bright neighbor
        # segment sits 2 pixels away, well inside the 0.4 arcsec
        # (~3.6 pixel) aperture.
        data = np.zeros(self.shape)
        segm = np.zeros(self.shape, dtype=np.int32)
        segm[40, 40] = 1
        segm[39:42, 42:44] = 2
        data[segm == 2] = 100.0
        model = SimpleNamespace(data=data << u.nJy, err=np.ones(self.shape) << u.nJy)
        xypos = np.array([[40.0, 40.0]])

        cat = ApertureCatalog(
            model, xypos, self._uniform_area_map(), SegmentationImage(segm), [1]
        )
        assert u.allclose(cat.aper04_flux, 0.0 * u.nJy, atol=1e-6 * u.nJy)

        # Without masking, the neighbor contributes to the aperture sum
        radius = cat.aperture_radii["circle_pix"][2, 0]
        unmasked = aperture_photometry(data, CircularAperture(xypos, radius))
        assert unmasked["aperture_sum"][0] > 100.0

    def test_neighbor_pixels_excluded_from_annulus_background(self):
        # Fill the right half of the annulus with a bright neighbor
        # segment. Its pixels dominate an unmasked median but must not
        # affect the masked one.
        data = np.zeros(self.shape)
        segm = np.zeros(self.shape, dtype=np.int32)
        segm[40, 40] = 1
        yy, xx = np.mgrid[: self.shape[0], : self.shape[1]]
        radius = np.hypot(xx - 40, yy - 40)
        r_in, r_out = np.array(
            ApertureCatalog.ANNULUS_RADII_ARCSEC
        ) / self.pixel_scale.to_value(u.arcsec)
        neighbor = (radius >= r_in - 1) & (radius <= r_out + 1) & (xx > 40)
        segm[neighbor] = 2
        data[neighbor] = 100.0
        model = SimpleNamespace(data=data << u.nJy, err=np.ones(self.shape) << u.nJy)
        xypos = np.array([[40.0, 40.0]])

        cat = ApertureCatalog(
            model, xypos, self._uniform_area_map(), SegmentationImage(segm), [1]
        )
        assert u.allclose(cat.aper_bkg_flux, 0.0 * cat.aper_bkg_flux.unit)

    def test_target_pixels_excluded_from_annulus_background(self):
        # The target's own segment extends through the right half of the
        # annulus with bright pixels. The local background must exclude
        # them, not just the pixels of other sources.
        data = np.zeros(self.shape)
        segm = np.zeros(self.shape, dtype=np.int32)
        yy, xx = np.mgrid[: self.shape[0], : self.shape[1]]
        radius = np.hypot(xx - 40, yy - 40)
        r_out = ApertureCatalog.ANNULUS_RADII_ARCSEC[1] / self.pixel_scale.to_value(
            u.arcsec
        )
        target = (radius <= r_out + 1) & (xx >= 40)
        segm[target] = 1
        data[target] = 100.0
        model = SimpleNamespace(data=data << u.nJy, err=np.ones(self.shape) << u.nJy)
        xypos = np.array([[40.0, 40.0]])

        cat = ApertureCatalog(
            model, xypos, self._uniform_area_map(), SegmentationImage(segm), [1]
        )
        assert u.allclose(cat.aper_bkg_flux, 0.0 * cat.aper_bkg_flux.unit)

    def test_matches_segmentation_catalog_circular_photometry(self):
        # Two overlapping Gaussian sources whose segments touch. The
        # aperture catalog must reproduce the photutils SourceCatalog
        # circular photometry with aperture_mask_method="mask".
        yy, xx = np.mgrid[: self.shape[0], : self.shape[1]]
        data = 200.0 * np.exp(-0.5 * ((xx - 40) ** 2 + (yy - 40) ** 2) / 2.0**2)
        data += 300.0 * np.exp(-0.5 * ((xx - 46) ** 2 + (yy - 41) ** 2) / 2.5**2)
        err = np.ones(self.shape)
        segm = np.zeros(self.shape, dtype=np.int32)
        segm[(data > 5) & (xx <= 43)] = 1
        segm[(data > 5) & (xx > 43)] = 2
        segm_img = SegmentationImage(segm)

        segm_cat = SourceCatalog(
            data << u.nJy, segm_img, error=err << u.nJy, aperture_mask_method="mask"
        )
        xypos = np.transpose((segm_cat.x_centroid, segm_cat.y_centroid))
        model = SimpleNamespace(data=data << u.nJy, err=err << u.nJy)
        cat = ApertureCatalog(
            model, xypos, self._uniform_area_map(), segm_img, segm_cat.labels
        )

        # A uniform pixel area map gives a single radius bin, so the
        # aperture radius is exact
        assert len(cat._radius_bins) == 1
        for name, radius in zip(
            cat.aperture_flux_colnames,
            cat.aperture_radii["circle_pix"][:, 0],
            strict=True,
        ):
            flux, flux_err = segm_cat.circular_photometry(radius)
            assert u.allclose(getattr(cat, name), flux, rtol=1e-6)
            assert u.allclose(getattr(cat, f"{name}_err"), flux_err, rtol=1e-6)

        # The masking changes the answer for the overlapping sources
        unmasked_cat = SourceCatalog(
            data << u.nJy, segm_img, error=err << u.nJy, aperture_mask_method="none"
        )
        radius = cat.aperture_radii["circle_pix"][2, 0]
        unmasked_flux, _ = unmasked_cat.circular_photometry(radius)
        assert not u.allclose(cat.aper04_flux, unmasked_flux, rtol=1e-3)
