from pathlib import Path
from re import match

import astropy.units as u
import numpy as np
import pyarrow
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from astropy.time import Time
from numpy.testing import assert_equal
from photutils.segmentation import SegmentationImage
from roman_datamodels import datamodels as rdm
from roman_datamodels import stnode
from roman_datamodels.datamodels import (
    ForcedImageSourceCatalogModel,
    ImageModel,
    ImageSourceCatalogModel,
    MosaicModel,
    MosaicSegmentationMapModel,
    MosaicSourceCatalogModel,
    SegmentationMapModel,
)

from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.source_catalog_step import SourceCatalogStep

from .helpers import compare_model_and_parquet_metadata


def make_test_image():
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
    err = np.zeros_like(data) + noise_scale

    return data, err


@pytest.fixture
def mosaic_model():
    defaults = {
        "meta": {
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
    model.meta.ref_file = stnode.RefFile.create_fake_data()
    model.meta.filename = "none"
    model.meta.cal_step = stnode.L3CalStep.create_fake_data()
    model.cal_logs = stnode.CalLogs.create_fake_data()
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.weight = 1.0 / err
    return model


@pytest.fixture
def image_model():
    model = ImageModel.create_fake_data(shape=(101, 101))
    model.meta.filename = "none"
    model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    model.meta.cal_logs = stnode.CalLogs.create_fake_data()
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.meta.photometry.conversion_megajanskys = (0.3324 * u.MJy / u.sr).value
    return model


def test_forced_catalog(image_model, function_jail, ignore_parquet_metadata_paths):
    output_filename = "force_cat.parquet"
    step = SourceCatalogStep()
    _ = step.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=5,
        npixels=10,
        save_results=True,
        output_file="source_cat.asdf",
    )
    result_force = step.call(
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

    assert Path(output_filename).exists()
    catalog = Table.read(output_filename)
    has_forced_fields = False
    for field in catalog.dtype.names:
        if "forced_" in field:
            has_forced_fields = True
    assert has_forced_fields

    compare_model_and_parquet_metadata(
        image_model, output_filename, ignore_parquet_metadata_paths
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
    output_filename = "test_cat.parquet"
    step = SourceCatalogStep()
    result = step.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    if save_results:
        assert Path(output_filename).exists()
        compare_model_and_parquet_metadata(
            image_model, output_filename, ignore_parquet_metadata_paths
        )
        cat = Table.read(output_filename)
    else:
        # FIXME: test output_filename doesn't exists but due to
        # https://github.com/spacetelescope/romancal/issues/1960 it always will
        cat = result.source_catalog
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
    step = SourceCatalogStep()
    mosaic_model.meta.filename = "test_coadd.asdf"
    output_filename = "test_cat.parquet"

    # Create model and set some crucial meta required to
    # create the L3 PSF for flux determination.
    result = step.call(
        mosaic_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    if save_results:
        assert Path(output_filename).exists()
        cat = Table.read(output_filename)
        compare_model_and_parquet_metadata(
            mosaic_model, output_filename, ignore_parquet_metadata_paths
        )
    else:
        # FIXME: test output_filename doesn't exists but due to
        # https://github.com/spacetelescope/romancal/issues/1960 it always will
        cat = result.source_catalog
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


def test_background(mosaic_model, function_jail):
    """
    Test background fallback when Background2D fails.
    """
    step = SourceCatalogStep()
    result = step.call(
        mosaic_model,
        bkg_boxsize=1000,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=25,
        fit_psf=False,
    )

    cat = result.source_catalog

    assert isinstance(cat, Table)


def test_l2_input_model_unchanged(image_model, function_jail):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    original_data = image_model.data.copy()
    original_err = image_model.err.copy()

    step = SourceCatalogStep()
    step.call(
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


def test_l3_input_model_unchanged(mosaic_model, function_jail):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    original_data = mosaic_model.data.copy()
    original_err = mosaic_model.err.copy()

    step = SourceCatalogStep()
    step.call(
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
        step = SourceCatalogStep()
        result = step.call(model)
        cat = result.source_catalog
        assert isinstance(cat, Table)
        assert len(cat) == 0


def test_inputs(mosaic_model):
    data = np.ones((3, 3), dtype=int)
    data[1, 1] = 1
    segm = SegmentationImage(data)
    cdata = np.ones((3, 3))
    kernel_fwhm = 2.0
    with pytest.raises(ValueError):
        RomanSourceCatalog(np.ones((3, 3)), segm, cdata, kernel_fwhm, fit_psf=True)


def test_psf_photometry(function_jail, image_model):
    """
    Test PSF photometry.
    """
    step = SourceCatalogStep()
    result = step.call(
        image_model,
        bkg_boxsize=20,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
        save_results=False,
    )

    cat = result.source_catalog
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
    step = SourceCatalogStep()
    result = step.call(
        image_model,
        bkg_boxsize=20,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
        save_results=False,
        fit_psf=fit_psf,
    )

    cat = result.source_catalog
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
    "snr_threshold, npixels, nsources, save_results, return_updated_model, expected_result, expected_outputs",
    (
        (
            3,
            10,
            7,
            True,
            True,
            ImageModel,
            {
                "cat": ImageSourceCatalogModel,
                "segm": SegmentationMapModel,
                "sourcecatalog": ImageModel,
            },
        ),
        (
            3,
            50,
            5,
            True,
            False,
            ImageSourceCatalogModel,
            {
                "cat": ImageSourceCatalogModel,
                "segm": SegmentationMapModel,
            },
        ),
        (
            10,
            10,
            7,
            False,
            True,
            ImageModel,
            {
                "cat": ImageSourceCatalogModel,
                "segm": SegmentationMapModel,
            },
        ),
        (
            20,
            10,
            5,
            False,
            False,
            ImageSourceCatalogModel,
            {
                "cat": ImageSourceCatalogModel,
                "segm": SegmentationMapModel,
            },
        ),
    ),
)
def test_l2_source_catalog_keywords(
    image_model,
    snr_threshold,
    npixels,
    nsources,
    save_results,
    return_updated_model,
    expected_result,
    expected_outputs,
    monkeypatch,
    function_jail,
):
    """
    Test that the proper object is returned in the call to SourceCatalogStep
    and that the desired output files are saved to the disk with the correct type.
    """
    monkeypatch.setattr(
        SourceCatalogStep, "return_updated_model", return_updated_model, raising=False
    )

    result = SourceCatalogStep.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    # assert that we returned the correct object
    assert isinstance(result, expected_result)

    # assert that the desired output files were saved to disk and that
    # they are of the correct type
    for suffix in expected_outputs.keys():
        if suffix == "cat":
            ext = "parquet"
        else:
            ext = "asdf"

        # annoying case.  Sometimes we have meta.filename as just "none" and
        # this test relies on the filename actually being at none_cat.parquet, etc.
        # But if we return a source catalog with a correct meta.filename (e.g.,
        # none_cat.parquet), this test needs to know how to translate that back
        # to the equivalent segmentation file.
        basefilename = result.meta.filename.split("_")[0]
        filepath = Path(function_jail / f"{basefilename}_{suffix}.{ext}")
        assert filepath.exists()

        if suffix == "cat":
            # the catalog is saved as a parquet file
            tbl = pyarrow.parquet.read_table(filepath)
            assert isinstance(tbl, pyarrow.Table)
        else:
            assert isinstance(rdm.open(filepath), expected_outputs.get(suffix))


@pytest.mark.parametrize(
    "snr_threshold, npixels, nsources, save_results, return_updated_model, expected_result, expected_outputs",
    (
        (
            3,
            10,
            7,
            True,
            True,
            MosaicModel,
            {
                "cat": MosaicSourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
                "sourcecatalog": MosaicModel,
            },
        ),
        (
            3,
            50,
            5,
            True,
            False,
            MosaicSourceCatalogModel,
            {
                "cat": MosaicSourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
            },
        ),
        (
            10,
            10,
            7,
            False,
            True,
            MosaicModel,
            {
                "cat": MosaicSourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
            },
        ),
        (
            20,
            10,
            5,
            False,
            False,
            MosaicSourceCatalogModel,
            {
                "cat": MosaicSourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
            },
        ),
    ),
)
def test_l3_source_catalog_keywords(
    mosaic_model,
    snr_threshold,
    npixels,
    nsources,
    save_results,
    return_updated_model,
    expected_result,
    expected_outputs,
    monkeypatch,
    function_jail,
):
    """
    Test that the proper object is returned in the call to SourceCatalogStep
    and that the desired output files are saved to the disk with the correct type.
    """
    # this step attribute controls whether to return a datamodel or source catalog
    monkeypatch.setattr(
        SourceCatalogStep, "return_updated_model", return_updated_model, raising=False
    )

    result = SourceCatalogStep.call(
        mosaic_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    # assert that we returned the correct object
    assert isinstance(result, expected_result)

    # assert that the desired output files were saved to disk and that
    # they are of the correct type
    for suffix in expected_outputs.keys():
        if suffix == "cat":
            ext = "parquet"
        else:
            ext = "asdf"

        basefilename = result.meta.filename.split("_")[0]
        filepath = Path(function_jail / f"{basefilename}_{suffix}.{ext}")
        assert filepath.exists()

        if suffix == "cat":
            # the catalog is saved as a parquet file
            tbl = pyarrow.parquet.read_table(filepath)
            assert isinstance(tbl, pyarrow.Table)
        else:
            assert isinstance(rdm.open(filepath), expected_outputs.get(suffix))


@pytest.mark.parametrize(
    "return_updated_model, expected_result",
    (
        (
            True,
            ImageModel,
        ),
        (
            False,
            ImageSourceCatalogModel,
        ),
    ),
)
def test_l2_source_catalog_return_updated_model_attribute(
    image_model,
    return_updated_model,
    expected_result,
    function_jail,
):
    """
    Test that the proper object is returned in the call to SourceCatalogStep.
    """
    step = SourceCatalogStep(
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=3,
        npixels=10,
    )

    if return_updated_model:
        # mimic what happens in the ELP -- i.e. set the "hidden" parameter
        # to cause this step to return a model instead of a catalog
        step.return_updated_model = return_updated_model

    result = step.run(image_model)

    # assert that we returned the correct object
    assert isinstance(result, expected_result)
