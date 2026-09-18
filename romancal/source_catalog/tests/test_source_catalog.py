from pathlib import Path
from re import match
from types import SimpleNamespace

import asdf
import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from astropy.time import Time
from numpy.testing import assert_allclose, assert_equal
from roman_datamodels import datamodels as rdm
from roman_datamodels.datamodels import (
    ForcedImageSourceCatalogModel,
    ImageModel,
    ImageSourceCatalogModel,
    MosaicModel,
    MosaicSegmentationMapModel,
    MosaicSourceCatalogModel,
    SegmentationMapModel,
)

from romancal.source_catalog._skyvals import compute_skyvals
from romancal.source_catalog._source_catalog import RomanSourceCatalog
from romancal.source_catalog._unit_conversion import (
    validate_and_convert_to_flux_density,
)
from romancal.source_catalog._wcs_utils import pixel_area_map
from romancal.source_catalog.source_catalog_step import SourceCatalogStep

from .helpers import compare_model_and_parquet_metadata

SKYVALS_DTYPE = np.dtype(
    [
        ("healpix17", np.int64),
        ("data", np.float32),
        ("err", np.float32),
        ("covfrac", np.float32),
    ]
)


def make_test_image(err_dtype=np.float16):
    g1 = Gaussian2D(121.0, 11.1, 12.2, 1.5, 1.5)
    g2 = Gaussian2D(70, 65, 18, 9.2, 4.5)
    g3 = Gaussian2D(111.0, 41, 42.7, 8.0, 3.0, theta=30 * u.deg)
    g4 = Gaussian2D(81.0, 17, 52.7, 4, 2, theta=102 * u.deg)
    g5 = Gaussian2D(107.0, 65, 71, 12, 2, theta=142 * u.deg)
    g6 = Gaussian2D(50, 20, 80, 2.1, 2.1)
    g7 = Gaussian2D(97.0, 85, 87.3, 4, 2, theta=-30 * u.deg)

    yy, xx = np.mgrid[0:101, 0:101]
    data = (
        g1(xx, yy)
        + g2(xx, yy)
        + g3(xx, yy)
        + g4(xx, yy)
        + g5(xx, yy)
        + g6(xx, yy)
        + g7(xx, yy)
    ).value.astype("float32")

    y0 = 2
    x0 = 90
    dd = 5
    value = -20
    y1 = y0 + dd
    x1 = x0 + dd - 1
    x2 = x0 + dd
    x3 = x0 + 2 * dd - 1
    data[y0:y1, x0] = value
    data[y0, x0:x1] = value
    data[y0:y1, x2] = value
    data[y0, x2:x3] = value
    data[y0 + dd // 2, x2:x3] = value
    data[y1 - 1, x2:x3] = value
    data[y0 + 1 : y1 - 1, x3] = value

    rng = np.random.default_rng(seed=123)
    noise_scale = 2.5
    noise = rng.normal(0, noise_scale, size=data.shape)
    data += noise
    err = (np.zeros_like(data) + noise_scale).astype(err_dtype)

    return data, err


@pytest.fixture
def mosaic_model():
    defaults = {
        "meta": {
            "data_release_id": "r1",
            "coadd_info": {
                "time_first": Time("2024-01-01T12:00:00.000", format="isot"),
            },
            "instrument": {
                "optical_element": "F158",
            },
            "resample": {"pixfrac": 1.0},
            "wcsinfo": {
                "pixel_scale": 1.5277777769528157e-05
            },  # Taken from regtest test L3 mosaic.
        }
    }
    model = MosaicModel.create_fake_data(defaults=defaults, shape=(101, 101))
    model.meta.filename = "none"
    model.meta.cal_step = {}
    for step_name in model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        model.meta.cal_step[step_name] = "INCOMPLETE"
    model.cal_logs = []
    data, err = make_test_image(err_dtype=np.float32)
    model.data = data
    model.err = err
    model.weight = 1.0 / err
    return model


@pytest.fixture
def image_model():
    model = ImageModel.create_fake_data(shape=(101, 101))
    model.meta.exposure.start_time = Time(
        "2024-01-03T00:00:00.0", format="isot", scale="utc"
    )
    model.meta.filename = "none"
    model.meta.cal_step = {}
    for step_name in model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        model.meta.cal_step[step_name] = "INCOMPLETE"
    model.meta.cal_logs = []
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.meta.photometry.conversion_megajanskys = (0.3324 * u.MJy / u.sr).value
    return model


def test_forced_catalog(image_model, function_jail, ignore_parquet_metadata_paths):
    """Purpose: forced photometry succeeds when the forcing segm has detection_image."""
    output_filename = "force_cat.parquet"
    _ = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=True,
        output_file="source_cat.asdf",
    )
    result_force, segmentation_map = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=True,
        output_file=output_filename,
        forced_segmentation="source_segm.asdf",
    )
    assert isinstance(result_force, ForcedImageSourceCatalogModel)
    assert isinstance(segmentation_map, SegmentationMapModel)

    assert Path(output_filename).exists()
    assert Path("force_segm.asdf").exists()
    catalog = Table.read(output_filename)
    has_forced_fields = False
    for field in catalog.dtype.names:
        if "forced_" in field:
            has_forced_fields = True
    assert has_forced_fields

    # The unprefixed columns are measured from the detection image saved
    # in the forcing segmentation file. The same image is used here, so
    # they must match the original detection catalog. The centroid
    # errors depend on the detection image flux scale.
    detection_catalog = Table.read("source_cat.parquet")
    assert_equal(catalog["label"], detection_catalog["label"])
    err_names = [
        name
        for name in catalog.colnames
        if "centroid" in name
        and name.endswith("_err")
        and not name.startswith("forced_")
    ]
    assert "x_centroid_err" in err_names
    assert "y_centroid_win_err" in err_names
    for name in err_names:
        assert_allclose(catalog[name], detection_catalog[name], rtol=1e-5)

    # The saved detection images are in flux density units and carry
    # a unit marker so that forced photometry does not convert them again
    with (
        rdm.open("source_segm.asdf") as source_segm,
        rdm.open("force_segm.asdf") as force_segm,
    ):
        assert source_segm.detection_image_unit == "nJy"
        assert force_segm.detection_image_unit == "nJy"
        expected = np.asarray(source_segm.detection_image)
        assert_equal(np.asarray(force_segm.detection_image), expected)
        assert_equal(np.asarray(segmentation_map.detection_image), expected)

    for name in ("x_centroid_err", "y_centroid_win_err", "ra_centroid_err"):
        assert np.all(np.isfinite(catalog[name]))
        assert np.all(catalog[name] > 0)

    compare_model_and_parquet_metadata(
        image_model, output_filename, ignore_parquet_metadata_paths
    )


class TestConvolvedDataUnits:
    """
    Test the conversion of convolved data whose units differ from the
    model data.
    """

    l2_to_sb = 2.0
    sb_to_flux = np.full((101, 101), 10.0) * u.nJy

    def convert(self, model, convolved_data):
        return validate_and_convert_to_flux_density(
            model,
            convolved_data,
            flux_unit=u.nJy,
            l2_to_sb=self.l2_to_sb,
            sb_to_flux=self.sb_to_flux,
        )

    def test_convolved_with_unit(self, image_model):
        """
        Convolved data already in flux density units must not be scaled
        with a unitless model.
        """
        data = image_model.data.copy()
        convolved_data = np.full((101, 101), 5.0, dtype=np.float32) << u.uJy
        result = self.convert(image_model, convolved_data)
        assert_allclose(result, 5000.0 * u.nJy)
        assert_allclose(image_model.data, data * 20.0 * u.nJy, rtol=1e-6)

    def test_model_with_unit(self, image_model):
        """
        Unitless convolved data must be converted from Level-2 units
        when the model is already in flux density units.
        """
        image_model["data"] = image_model.data << u.uJy
        image_model["err"] = image_model.err.astype(np.float32) << u.uJy
        data = image_model.data.copy()
        convolved_data = np.full((101, 101), 5.0, dtype=np.float32)
        result = self.convert(image_model, convolved_data)
        assert_allclose(result, 100.0 * u.nJy)
        assert_allclose(image_model.data, data)
        assert image_model.data.unit == u.nJy


def _write_forcing_segm(image_model, filename, *, unit, scale=1.0):
    """
    Write a forcing segmentation file with a modified detection image
    unit key. A `None` unit removes the key.
    """
    SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=True,
        output_file="source_cat.asdf",
    )
    with asdf.open("source_segm.asdf", memmap=False, lazy_load=False) as af:
        roman = af.tree["roman"]
        roman["detection_image"] = np.asarray(roman["detection_image"]) * scale
        if unit is None:
            del roman["detection_image_unit"]
        else:
            roman["detection_image_unit"] = unit
        af.write_to(filename)


def _call_forced(image_model, filename):
    return SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=False,
        forced_segmentation=filename,
    )


@pytest.mark.parametrize(("unit", "scale"), [("uJy", 1e-3), (None, 1.0)])
def test_forced_catalog_detection_image_unit(
    image_model, function_jail, caplog, unit, scale
):
    """
    Test a forcing detection image in an equivalent unit and one from a
    file without the unit key, which is assumed to be in nJy.
    """
    filename = "modified_segm.asdf"
    _write_forcing_segm(image_model, filename, unit=unit, scale=scale)
    forced_cat, forced_segm = _call_forced(image_model, filename)

    assert ("Assuming its detection_image is in nJy" in caplog.text) == (unit is None)

    detection_catalog = Table.read("source_cat.parquet")
    for name in ("x_centroid_err", "y_centroid_win_err"):
        assert_allclose(
            forced_cat.source_catalog[name], detection_catalog[name], rtol=1e-5
        )

    assert forced_segm.detection_image_unit == "nJy"
    with rdm.open("source_segm.asdf") as source_segm:
        assert_allclose(
            forced_segm.detection_image,
            np.asarray(source_segm.detection_image),
            rtol=1e-6,
        )


@pytest.mark.parametrize(
    ("unit", "match"),
    [("s", "not equivalent to the desired flux unit"), ("bad", "not a valid unit")],
)
def test_forced_catalog_invalid_detection_image_unit(
    image_model, function_jail, unit, match
):
    filename = "modified_segm.asdf"
    _write_forcing_segm(image_model, filename, unit=unit)
    with pytest.raises(ValueError, match=match):
        _call_forced(image_model, filename)


def test_forced_catalog_requires_detection_image(image_model, function_jail):
    """Purpose: forced photometry errors clearly if forcing segm lacks detection_image."""
    _, segm = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=False,
    )
    # Simulate legacy / empty products that only carry the label map.
    bare_segm = SegmentationMapModel.create_minimal({"meta": segm.meta})
    bare_segm.data = segm.data.copy()
    forced_segm_path = Path("no_detection_segm.asdf")
    bare_segm.save(forced_segm_path)

    with pytest.raises(ValueError, match="must include a detection_image"):
        SourceCatalogStep.call(
            image_model,
            bkg_boxsize=50,
            kernel_fwhm=2.0,
            snr_threshold=5,
            npixels=10,
            save_results=False,
            forced_segmentation=str(forced_segm_path),
        )


@pytest.mark.parametrize(
    "snr_threshold, npixels, nsources, save_results",
    (
        (3, 10, 7, True),
        (3, 50, 5, False),
        (10, 10, 7, False),
        (20, 10, 5, False),
        (25, 10, 3, False),
        (35, 10, 1, False),
        (50, 10, 0, False),
    ),
)
def test_l2_source_catalog(
    image_model,
    snr_threshold,
    npixels,
    nsources,
    save_results,
    function_jail,
    ignore_parquet_metadata_paths,
):
    image_model.meta.filename = "test_cal.asdf"
    catalog_filename = "test_cat.parquet"
    segmentation_map_filename = "test_segm.asdf"

    result_catalog, result_segmentation_map = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    assert isinstance(result_catalog, ImageSourceCatalogModel)
    assert isinstance(result_segmentation_map, SegmentationMapModel)

    if save_results:
        assert Path(segmentation_map_filename).exists()
        assert Path(catalog_filename).exists()
        compare_model_and_parquet_metadata(
            image_model, catalog_filename, ignore_parquet_metadata_paths
        )
        cat = Table.read(catalog_filename)
    else:
        assert not Path(segmentation_map_filename).exists()
        assert not Path(catalog_filename).exists()
        cat = result_catalog.source_catalog
        assert isinstance(cat, Table)
    assert len(cat) == nsources

    # Check that the ee_fraction_xx entries are in the metadata
    if "aperture_radii" in cat.meta:
        assert len(cat.meta["aperture_radii"]["circle_pix"]) > 0
        assert sum(1 for name in cat.colnames if match(r"^aper\d+_flux$", name)) == len(
            cat.meta["aperture_radii"]["circle_pix"]
        )

        assert "ee_fractions" in cat.meta
        assert len(cat.meta["ee_fractions"]) == len(
            cat.meta["aperture_radii"]["circle_pix"]
        )
    else:
        assert nsources == 0

    if nsources > 0:
        for colname in cat.colnames:
            if (
                "flux" in colname
                and "fluxfrac" not in colname
                and "aper_bkg_flux" not in colname
            ):
                assert cat[colname].unit == "nJy"
        assert np.min(cat["x_centroid"]) > 0.0
        assert np.min(cat["y_centroid"]) > 0.0
        assert np.max(cat["x_centroid"]) < 100.0
        assert np.max(cat["y_centroid"]) < 100.0
        assert np.any(cat["ra"])
        assert np.any(cat["dec"])


@pytest.mark.parametrize(
    "snr_threshold, npixels, nsources, save_results",
    (
        (3, 10, 7, True),
        (3, 50, 5, False),
        (10, 10, 7, False),
        (20, 10, 5, False),
        (25, 10, 3, False),
        (35, 10, 1, False),
        (50, 10, 0, False),
    ),
)
def test_l3_source_catalog(
    mosaic_model,
    snr_threshold,
    npixels,
    nsources,
    save_results,
    function_jail,
    ignore_parquet_metadata_paths,
):
    mosaic_model.meta.filename = "test_coadd.asdf"
    catalog_filename = "test_cat.parquet"
    segmentation_map_filename = "test_segm.asdf"

    # Create model and set some crucial meta required to
    # create the L3 PSF for flux determination.
    result_catalog, result_segmentation_map = SourceCatalogStep.call(
        mosaic_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    assert isinstance(result_catalog, MosaicSourceCatalogModel)
    assert isinstance(result_segmentation_map, MosaicSegmentationMapModel)

    if save_results:
        assert Path(segmentation_map_filename).exists()
        assert Path(catalog_filename).exists()
        cat = Table.read(catalog_filename)
        compare_model_and_parquet_metadata(
            mosaic_model, catalog_filename, ignore_parquet_metadata_paths
        )
    else:
        assert not Path(segmentation_map_filename).exists()
        assert not Path(catalog_filename).exists()
        cat = result_catalog.source_catalog
        assert isinstance(cat, Table)
    assert len(cat) == nsources

    assert result_catalog.meta.data_release_id == mosaic_model.meta.data_release_id

    # Check that the ee_fraction_xx entries are in the metadata
    if "aperture_radii" in cat.meta:
        assert len(cat.meta["aperture_radii"]["circle_pix"]) > 0
        assert sum(1 for name in cat.colnames if match(r"^aper\d+_flux$", name)) == len(
            cat.meta["aperture_radii"]["circle_pix"]
        )

        assert "ee_fractions" in cat.meta
        assert len(cat.meta["ee_fractions"]) == len(
            cat.meta["aperture_radii"]["circle_pix"]
        )
    else:
        assert nsources == 0

    if nsources > 0:
        for colname in cat.colnames:
            if (
                "flux" in colname
                and "fluxfrac" not in colname
                and "aper_bkg_flux" not in colname
            ):
                assert cat[colname].unit == "nJy"
        assert np.min(cat["x_centroid"]) > 0.0
        assert np.min(cat["y_centroid"]) > 0.0
        assert np.max(cat["x_centroid"]) < 100.0
        assert np.max(cat["y_centroid"]) < 100.0
        assert np.any(cat["ra"])
        assert np.any(cat["dec"])


def test_centroid_errors_and_sky_orientation(image_model, function_jail):
    """
    The centroid errors, sky orientation, and annulus background error
    are calculated by photutils and must be populated with valid values.
    """
    result_catalog, _ = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
        save_results=False,
    )
    cat = result_catalog.source_catalog
    assert len(cat) > 0

    pix_names = (
        "x_centroid_err",
        "y_centroid_err",
        "x_centroid_win_err",
        "y_centroid_win_err",
    )
    sky_names = (
        "ra_centroid_err",
        "dec_centroid_err",
        "ra_centroid_win_err",
        "dec_centroid_win_err",
    )
    for names, unit in ((pix_names, u.pix), (sky_names, u.arcsec)):
        for name in names:
            assert cat[name].unit == unit
            assert cat[name].dtype == np.float32
            assert np.all(np.isfinite(cat[name]))
            assert np.all(cat[name] > 0)

    # The total sky error is the total pixel error times the pixel
    # scale, independent of the WCS rotation
    wcs = image_model.meta.wcs
    pixel_scale = np.sqrt(pixel_area_map(wcs, image_model.data.shape).mean())
    pix_err = np.hypot(cat["x_centroid_err"].value, cat["y_centroid_err"].value)
    sky_err = np.hypot(cat["ra_centroid_err"].value, cat["dec_centroid_err"].value)
    assert_allclose(sky_err / pix_err, pixel_scale.to_value(u.arcsec), rtol=1e-3)

    assert cat["orientation_sky"].unit == u.deg
    assert cat["orientation_sky"].dtype == np.float32
    assert np.all(cat["orientation_sky"] > -90)
    assert np.all(cat["orientation_sky"] <= 90)

    assert cat["aper_bkg_flux_err"].dtype == np.float32
    assert np.all(cat["aper_bkg_flux_err"] > 0)


def test_background(mosaic_model, function_jail):
    """
    Test background fallback when Background2D fails.
    """
    result_catalog, _ = SourceCatalogStep.call(
        mosaic_model,
        bkg_boxsize=1000,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=25,
        fit_psf=False,
    )

    cat = result_catalog.source_catalog

    assert isinstance(cat, Table)


@pytest.mark.parametrize("model_fixture", ("image_model", "mosaic_model"))
def test_source_catalog_populates_dust_ebv(model_fixture, request, function_jail):
    """Ensure prompt source catalogs include a per-source dust_ebv column."""
    model = request.getfixturevalue(model_fixture)
    result_catalog, _ = SourceCatalogStep.call(
        model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
        save_results=False,
    )
    cat = result_catalog.source_catalog
    assert "dust_ebv" in cat.colnames
    assert len(cat["dust_ebv"]) == len(cat)
    assert cat["dust_ebv"].dtype == np.float32


def test_nested_metadata_propagated_to_catalog_and_segmentation(
    image_model, function_jail
):
    """
    Nested (list-of-list) metadata such as ``meta.exposure.read_pattern``
    should be propagated to both the source catalog and segmentation map
    models.
    """
    read_pattern = [[1], [2, 3], [4, 5, 6], [7, 8, 9, 10]]
    image_model.meta.exposure.read_pattern = read_pattern

    result_catalog, result_segmentation_map = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=False,
        fit_psf=False,
    )

    assert isinstance(result_catalog, ImageSourceCatalogModel)
    assert isinstance(result_segmentation_map, SegmentationMapModel)

    for result in (result_catalog, result_segmentation_map):
        assert [list(r) for r in result.meta.exposure.read_pattern] == read_pattern


def test_l2_input_model_unchanged(image_model, function_jail):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    original_data = image_model.data.copy()
    original_err = image_model.err.copy()

    SourceCatalogStep.call(
        image_model,
        snr_threshold=0.5,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
        fit_psf=False,
    )

    assert_equal(original_data, image_model.data)
    assert_equal(original_err, image_model.err)


def test_l2_segmentation_contains_skyvals(image_model):
    _, result_segmentation_map = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=False,
        fit_psf=False,
    )

    assert isinstance(result_segmentation_map, SegmentationMapModel)
    assert "skyvals" in result_segmentation_map
    assert "healpix11_cov" in result_segmentation_map

    skyvals = result_segmentation_map.skyvals
    assert skyvals.dtype == SKYVALS_DTYPE
    assert skyvals.shape[0] > 0

    healpix11_cov = result_segmentation_map.healpix11_cov
    assert healpix11_cov.ndim == 1
    assert healpix11_cov.dtype == np.int64


def test_l2_segmentation_without_skyvals_when_disabled(image_model):
    _, result_segmentation_map = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=False,
        fit_psf=False,
        compute_skyvals=False,
    )

    assert isinstance(result_segmentation_map, SegmentationMapModel)
    assert "skyvals" not in result_segmentation_map
    assert "healpix11_cov" not in result_segmentation_map


def test_skyvals_coverage_includes_source_pixels(monkeypatch):
    class WCS:
        def pixel_to_world_values(self, x, y):
            return np.asarray(x) * 180.0, np.asarray(y) * 0.0

    input_model = SimpleNamespace(
        data=np.array([[10.0, 20.0]], dtype=np.float32),
        err=np.ones((1, 2), dtype=np.float32),
        meta=SimpleNamespace(wcs=WCS()),
    )
    segmentation = np.array([[1, 0]], dtype=np.uint32)
    bad_pixel_mask = np.zeros((1, 2), dtype=bool)

    monkeypatch.setattr(
        "romancal.source_catalog._skyvals.get_pixel_area_sr", lambda model: 1.0
    )
    skyvals, healpix11_cov = compute_skyvals(input_model, segmentation, bad_pixel_mask)

    assert skyvals.shape == (1,)
    assert skyvals["data"][0] == 20.0
    assert healpix11_cov.shape == (2,)


def test_l2_skyvals_values_and_covfrac_reasonable(image_model):
    """
    Verify skyvals statistics are sensible on a controlled input image.
    """
    rng = np.random.default_rng(seed=9)
    image_model.data = (rng.normal(0, 1, size=image_model.data.shape) + 100.0).astype(
        np.float32
    )
    image_model.err = np.full_like(image_model.data, 3.0, dtype=np.float32)

    _, result_segmentation_map = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=False,
        fit_psf=False,
    )

    skyvals = result_segmentation_map.skyvals
    covfrac = skyvals["covfrac"]
    assert skyvals.shape[0] > 0
    # data medians should recover the injected background offset
    # (atol=5.0 allows for some variation due to the random noise)
    np.testing.assert_allclose(np.nanmedian(skyvals["data"]), 100.0, atol=5.0)
    # covfrac should be bounded in [0, 1]
    assert np.all(covfrac >= 0.0)
    assert np.all(covfrac <= 1.0)
    # at least some healpixels should have near-full coverage
    assert np.any(covfrac > 0.9)


def test_l3_input_model_unchanged(mosaic_model, function_jail):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    original_data = mosaic_model.data.copy()
    original_err = mosaic_model.err.copy()

    SourceCatalogStep.call(
        mosaic_model,
        snr_threshold=0.5,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        save_results=False,
        fit_psf=False,
    )

    assert_equal(original_data, mosaic_model.data)
    assert_equal(original_err, mosaic_model.err)


def test_invalid_step_inputs(image_model, mosaic_model, function_jail):
    for input_model in (image_model, mosaic_model):
        model = input_model.copy()
        model.data = np.full(model.data.shape, np.nan)
        result_catalog, _ = SourceCatalogStep.call(model)
        cat = result_catalog.source_catalog
        assert isinstance(cat, Table)
        assert len(cat) == 0


def test_inputs(mosaic_model):
    with pytest.raises(ValueError, match="The input model must be an"):
        RomanSourceCatalog(None, None, None, None, None)


def test_psf_photometry(function_jail, image_model):
    """
    Test PSF photometry.
    """
    result_catalog, _ = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=20,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
        save_results=False,
    )

    cat = result_catalog.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7

    for colname in cat.colnames:
        if (
            "flux" in colname
            and "fluxfrac" not in colname
            and "aper_bkg_flux" not in colname
        ):
            assert cat[colname].unit == "nJy"

    for colname in cat.colnames:
        if "psf" in colname:
            assert len(cat[colname])  # make sure the column isn't empty
            assert not np.any(np.isnan(cat[colname]))  # and contains no nans


@pytest.mark.parametrize("fit_psf", [True, False])
def test_do_psf_photometry_column_names(function_jail, image_model, fit_psf):
    """
    Test that fit_psf will determine whether the PSF
    photometry columns are added to the final catalog or not.
    """
    result_catalog, _ = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=20,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
        save_results=False,
        fit_psf=fit_psf,
    )

    cat = result_catalog.source_catalog
    assert isinstance(cat, Table)

    psf_colnames = []
    for colname in cat.colnames:
        if "psf" in colname:
            psf_colnames.append(colname)

    if fit_psf:
        assert len(psf_colnames) > 0
    else:
        assert len(psf_colnames) == 0


@pytest.mark.parametrize(
    "ra, dec",
    [
        (np.array([10.0, 20.0]), np.array([30.0])),
        (np.array([[10.0, 20.0]]), np.array([30.0, 40.0])),
    ],
)
def test_get_dust_ebv_shape_mismatch_raises(ra, dec):
    """Raise when RA/Dec input shapes are inconsistent."""
    cat = object.__new__(RomanSourceCatalog)
    map_paths = {
        RomanSourceCatalog.north_galactic_pole_id: "north.fits",
        RomanSourceCatalog.south_galactic_pole_id: "south.fits",
    }
    with pytest.raises(ValueError, match=r"ra\.shape must equal dec\.shape"):
        cat._get_dust_ebv(ra, dec, map_paths)


def test_dust_ebv_property_returns_nan_on_failure(monkeypatch):
    """Return NaNs when CRDS lookup/interpolation fails."""

    def fail_getreferences(*args, **kwargs):
        raise RuntimeError()

    monkeypatch.setattr(
        "romancal.source_catalog._source_catalog.getreferences", fail_getreferences
    )

    cat = object.__new__(RomanSourceCatalog)
    cat.ra = np.array([1.0, 2.0, 3.0], dtype=float)
    cat.dec = np.array([4.0, 5.0, 6.0], dtype=float)
    cat.n_sources = 3

    result = cat.dust_ebv
    assert result.dtype == np.float32
    assert result.shape == (3,)
    assert np.all(np.isnan(result))
