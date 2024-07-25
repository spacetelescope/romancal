import os
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from numpy.testing import assert_allclose
from photutils.segmentation import SegmentationImage
from roman_datamodels import datamodels as rdm
from roman_datamodels.datamodels import (
    ImageModel,
    MosaicModel,
    MosaicSegmentationMapModel,
    MosaicSourceCatalogModel,
    SourceCatalogModel,
)
from roman_datamodels.maker_utils import mk_level2_image, mk_level3_mosaic

from romancal.source_catalog.reference_data import ReferenceData
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.source_catalog_step import SourceCatalogStep


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
    noise = rng.normal(0, 2.5, size=data.shape)
    data += noise
    err = data / 10.0

    return data, err


@pytest.fixture
def mosaic_model():
    wfi_mosaic = mk_level3_mosaic(shape=(101, 101))
    model = MosaicModel(wfi_mosaic)
    data, err = make_test_image()
    units = u.MJy / u.sr
    data <<= units
    err <<= units
    model.data = data
    model.err = err
    model.weight = 1.0 / err.value
    return model


@pytest.fixture
def image_model():
    wfi_image = mk_level2_image(shape=(101, 101))
    model = ImageModel(wfi_image)
    data, err = make_test_image()
    units = u.DN / u.s
    data <<= units
    err <<= units
    model.data = data
    model.err = err
    model.meta.photometry.conversion_megajanskys = 0.3324 * u.MJy / u.sr
    return model


@pytest.mark.webbpsf
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
    image_model, snr_threshold, npixels, nsources, save_results, tmp_path
):
    os.chdir(tmp_path)
    step = SourceCatalogStep()
    result = step.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    cat = result.source_catalog

    assert isinstance(cat, Table)

    columns = [
        "label",
        "xcentroid",
        "ycentroid",
        "ra_centroid",
        "dec_centroid",
        "aper_bkg_flux",
        "aper_bkg_flux_err",
        "aper30_flux",
        "aper30_flux_err",
        "aper50_flux",
        "aper50_flux_err",
        "aper70_flux",
        "aper70_flux_err",
        "aper_total_flux",
        "aper_total_flux_err",
        "CI_50_30",
        "CI_70_50",
        "CI_70_30",
        "is_extended",
        "sharpness",
        "roundness",
        "nn_label",
        "nn_dist",
        "isophotal_flux",
        "isophotal_flux_err",
        "isophotal_area",
        "semimajor_sigma",
        "semiminor_sigma",
        "ellipticity",
        "orientation",
        "sky_orientation",
        "ra_bbox_ll",
        "dec_bbox_ll",
        "ra_bbox_ul",
        "dec_bbox_ul",
        "ra_bbox_lr",
        "dec_bbox_lr",
        "ra_bbox_ur",
        "dec_bbox_ur",
    ]

    assert len(cat) == nsources

    if nsources > 0:
        for col in columns:
            assert col in cat.colnames
        assert np.min(cat["xcentroid"]) > 0.0
        assert np.min(cat["ycentroid"]) > 0.0
        assert np.max(cat["xcentroid"]) < 100.0
        assert np.max(cat["ycentroid"]) < 100.0


@pytest.mark.webbpsf
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
    mosaic_model, snr_threshold, npixels, nsources, save_results, tmp_path
):
    os.chdir(tmp_path)
    step = SourceCatalogStep()

    result = step.call(
        mosaic_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    cat = result.source_catalog

    assert isinstance(cat, Table)

    columns = [
        "label",
        "xcentroid",
        "ycentroid",
        "ra_centroid",
        "dec_centroid",
        "aper_bkg_flux",
        "aper_bkg_flux_err",
        "aper30_flux",
        "aper30_flux_err",
        "aper50_flux",
        "aper50_flux_err",
        "aper70_flux",
        "aper70_flux_err",
        "aper_total_flux",
        "aper_total_flux_err",
        "CI_50_30",
        "CI_70_50",
        "CI_70_30",
        "is_extended",
        "sharpness",
        "roundness",
        "nn_label",
        "nn_dist",
        "isophotal_flux",
        "isophotal_flux_err",
        "isophotal_area",
        "semimajor_sigma",
        "semiminor_sigma",
        "ellipticity",
        "orientation",
        "sky_orientation",
        "ra_bbox_ll",
        "dec_bbox_ll",
        "ra_bbox_ul",
        "dec_bbox_ul",
        "ra_bbox_lr",
        "dec_bbox_lr",
        "ra_bbox_ur",
        "dec_bbox_ur",
    ]

    assert len(cat) == nsources

    if nsources > 0:
        for col in columns:
            assert col in cat.colnames
        assert np.min(cat["xcentroid"]) > 0.0
        assert np.min(cat["ycentroid"]) > 0.0
        assert np.max(cat["xcentroid"]) < 100.0
        assert np.max(cat["ycentroid"]) < 100.0


@pytest.mark.webbpsf
def test_background(mosaic_model, tmp_path):
    """
    Test background fallback when Background2D fails.
    """
    os.chdir(tmp_path)
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


@pytest.mark.webbpsf
def test_l2_input_model_unchanged(image_model, tmp_path):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    os.chdir(tmp_path)
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
    )

    assert_allclose(original_data, image_model.data, atol=5.0e-5)
    assert_allclose(original_err, image_model.err, atol=5.0e-5)


@pytest.mark.webbpsf
def test_l3_input_model_unchanged(mosaic_model, tmp_path):
    """
    Test that the input model data and error arrays are unchanged after
    processing by SourceCatalogStep.
    """
    os.chdir(tmp_path)
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
    )

    assert_allclose(original_data, mosaic_model.data, atol=5.0e-5)
    assert_allclose(original_err, mosaic_model.err, atol=5.0e-5)


@pytest.mark.webbpsf
def test_inputs(mosaic_model):
    with pytest.raises(ValueError):
        ReferenceData(np.ones((3, 3)), (30, 50, 70))
    with pytest.raises(ValueError):
        aperture_ee = (70, 50, 30)
        ReferenceData(mosaic_model, aperture_ee)
    with pytest.raises(ValueError):
        aperture_ee = (30, 50)
        ReferenceData(mosaic_model, aperture_ee)
    with pytest.raises(ValueError):
        aperture_ee = (-1, 50, 70)
        ReferenceData(mosaic_model, aperture_ee)
    with pytest.raises(ValueError):
        aperture_ee = (40, 70, 150)
        ReferenceData(mosaic_model, aperture_ee)

    data = np.ones((3, 3), dtype=int)
    data[1, 1] = 1
    segm = SegmentationImage(data)
    cdata = np.ones((3, 3))
    aper_params = {}
    ci_thresh = 100.0
    with pytest.raises(ValueError):
        RomanSourceCatalog(
            np.ones((3, 3)), segm, cdata, aper_params, ci_thresh, 2.0, fit_psf=True
        )

    with pytest.raises(ValueError):
        RomanSourceCatalog(
            mosaic_model, segm, cdata, aper_params, (1.0, 2.0, 3.0), 2.0, fit_psf=True
        )


@pytest.mark.webbpsf
def test_do_psf_photometry(tmp_path, image_model):
    """
    Test that do_psf_photometry can recover mock sources and their position and photometry.
    """
    os.chdir(tmp_path)

    # get column names mapping for PSF photometry
    psf_colnames_mapping = (
        RomanSourceCatalog.get_psf_photometry_catalog_colnames_mapping()
    )
    psf_colnames = [
        x.get("new_name")
        for x in psf_colnames_mapping
        if x.get("old_name") in ["x_fit", "y_fit", "flux_fit"]
    ]

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

    # check the number of sources that have been detected
    assert len(cat) == 7
    # check that all sources have both position and flux determined (ignore errors/flags)
    for col_name in psf_colnames:
        assert len(cat[col_name])  # make sure the column isn't empty
        assert not np.any(np.isnan(cat[col_name]))  # and contains no nans


@pytest.mark.webbpsf
@pytest.mark.parametrize("fit_psf", [True, False])
def test_do_psf_photometry_column_names(tmp_path, image_model, fit_psf):
    """
    Test that fit_psf will determine whether the PSF
    photometry columns are added to the final catalog or not.
    """
    os.chdir(tmp_path)

    # get column names mapping for PSF photometry
    psf_colnames_mapping = (
        RomanSourceCatalog.get_psf_photometry_catalog_colnames_mapping()
    )

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

    # check if the PSF photometry column names are present or not based on fit_psf value
    psf_colnames_present = all(
        x.get("new_name") in cat.colnames for x in psf_colnames_mapping
    )
    psf_colnames_not_present = all(
        x.get("new_name") not in cat.colnames for x in psf_colnames_mapping
    )

    assert (fit_psf and psf_colnames_present) or (
        not fit_psf and psf_colnames_not_present
    )


@pytest.mark.webbpsf
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
                "cat": SourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
                "sourcecatalog": ImageModel,
            },
        ),
        (
            3,
            50,
            5,
            True,
            False,
            SourceCatalogModel,
            {
                "cat": SourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
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
                "cat": SourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
            },
        ),
        (
            20,
            10,
            5,
            False,
            False,
            SourceCatalogModel,
            {
                "cat": SourceCatalogModel,
                "segm": MosaicSegmentationMapModel,
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
    tmp_path,
):
    """
    Test that the proper object is returned in the call to SourceCatalogStep
    and that the desired output files are saved to the disk with the correct type.
    """
    os.chdir(tmp_path)
    step = SourceCatalogStep
    # this step attribute controls whether to return a datamodel or source catalog
    step.return_updated_model = return_updated_model

    result = step.call(
        image_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    # assert that we returned the correct object
    assert isinstance(result, expected_result)

    # assert that the desired output files were saved to disk
    assert all(
        Path(tmp_path / f"{result.meta.filename}_{suffix}.asdf").exists()
        for suffix in expected_outputs.keys()
    )

    # assert that the desired output files were saved with the correct datamodel type
    assert all(
        isinstance(
            rdm.open(Path(tmp_path / f"{result.meta.filename}_{suffix}.asdf")),
            expected_outputs.get(suffix),
        )
        for suffix in expected_outputs.keys()
    )


@pytest.mark.webbpsf
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
    tmp_path,
):
    """
    Test that the proper object is returned in the call to SourceCatalogStep
    and that the desired output files are saved to the disk with the correct type.
    """
    os.chdir(tmp_path)
    step = SourceCatalogStep
    # this step attribute controls whether to return a datamodel or source catalog
    step.return_updated_model = return_updated_model

    result = step.call(
        mosaic_model,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        snr_threshold=snr_threshold,
        npixels=npixels,
        save_results=save_results,
    )

    # assert that we returned the correct object
    assert isinstance(result, expected_result)

    # assert that the desired output files were saved to disk
    assert all(
        Path(tmp_path / f"{result.meta.filename}_{suffix}.asdf").exists()
        for suffix in expected_outputs.keys()
    )

    # assert that the desired output files were saved with the correct datamodel type
    assert all(
        isinstance(
            rdm.open(Path(tmp_path / f"{result.meta.filename}_{suffix}.asdf")),
            expected_outputs.get(suffix),
        )
        for suffix in expected_outputs.keys()
    )


@pytest.mark.webbpsf
@pytest.mark.parametrize(
    "return_updated_model, expected_result",
    (
        (
            True,
            ImageModel,
        ),
        (
            False,
            SourceCatalogModel,
        ),
    ),
)
def test_l2_source_catalog_return_updated_model_attribute(
    image_model,
    return_updated_model,
    expected_result,
    tmp_path,
):
    """
    Test that the proper object is returned in the call to SourceCatalogStep.
    """
    os.chdir(tmp_path)

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
