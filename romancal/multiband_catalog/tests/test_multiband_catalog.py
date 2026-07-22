from copy import deepcopy
from pathlib import Path
from re import match

import astropy.units as u
import numpy as np
import pyarrow
import pytest
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from astropy.time import Time
from roman_datamodels import datamodels as rdm
from roman_datamodels.datamodels import MosaicModel, MultibandSegmentationMapModel

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog import MultibandCatalogStep
from romancal.multiband_catalog._multiband_catalog import match_recovered_sources, make_source_grid
from romancal.skycell.tests.test_skycell_match import mk_gwcs

SI_SCALE = 5/4
RNG_SEED = 42
EXPTIME = 300
MEANFLUX = 0.2
RA_REF = 270.0 * u.deg
DEC_REF = 66.0 * u.deg
ROLL_REF = 0.0 * u.deg
SHAPE = (500,500)


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

    rng = np.random.default_rng(seed=42)
    noise_scale = 2.5
    noise = rng.normal(0, noise_scale, size=data.shape)
    data += noise
    err = np.zeros_like(data) + noise_scale

    return data, err


def make_si_test_image():
    flux_scale = MEANFLUX
    scale = SI_SCALE

    g1 = Gaussian2D(flux_scale * 121.0, scale * 11, scale * 12, 1.5, 1.5)
    g2 = Gaussian2D(flux_scale * 70, scale * 65, scale * 18, 9.2, 4.5)
    g3 = Gaussian2D(flux_scale * 111.0, scale * 41, scale * 43, 8.0, 3.0, theta=30 * u.deg)
    g4 = Gaussian2D(flux_scale * 81.0, scale * 17, scale * 53, 4, 2, theta=102 * u.deg)
    g5 = Gaussian2D(flux_scale * 107.0, scale * 65, scale * 71, 12, 2, theta=142 * u.deg)
    g6 = Gaussian2D(flux_scale * 50, scale * 20, scale * 80, 2.1, 2.1)
    g7 = Gaussian2D(flux_scale * 97.0, scale * 85, scale * 88, 4, 2, theta=-30 * u.deg)

    yy, xx = np.mgrid[0:125, 0:125]
    smalldata = np.zeros(shape=(len(yy), len(xx)))
    smalldata = (
        g1(xx, yy)
        + g2(xx, yy)
        + g3(xx, yy)
        + g4(xx, yy)
        + g5(xx, yy)
        + g6(xx, yy)
        + g7(xx, yy)
    ).value.astype("float32")

    # Rotate pattern to make the test image a little less regular
    data = np.tile(smalldata, (4, 1))
    data = np.append(data, np.tile(np.rot90(smalldata), (4, 1)), axis=1)
    data = np.append(data, np.tile(np.rot90(smalldata, k=2), (4, 1)), axis=1)
    data = np.append(data, np.tile(np.rot90(smalldata, k=-1), (4, 1)), axis=1)

    assert data.shape == (500, 500)

    rng = np.random.default_rng(seed=RNG_SEED)
    noise_scale = MEANFLUX / EXPTIME
    noise = rng.normal(loc=MEANFLUX, scale=noise_scale, size=data.shape)
    data += noise
    err =  np.sqrt(np.ones_like(data))

    return data, err


@pytest.fixture
def mosaic_model(shape=(101, 101)):
    model = MosaicModel.create_fake_data(shape=shape)
    data, err = make_test_image()
    model.data = data
    model.err = err
    model.var_rnoise = err**2
    model.weight = 1.0 / err
    model.meta.instrument.optical_element = "F184"
    model.meta.coadd_info.time_first = Time("2027-01-01T00:00:00")
    model.meta.wcsinfo.pixel_scale = 0.11 / 3600  # degrees

    model.meta.wcsinfo.ra_ref = 270.0  # degrees
    model.meta.wcsinfo.dec_ref = 66.0  # degrees
    model.meta.wcsinfo.roll_ref = 0.0  # degrees
    model.meta.coadd_info.exposure_time = EXPTIME  # seconds

    model.meta.resample.pixfrac = 0.5
    model.meta.data_release_id = "r1"
    return model


@pytest.fixture
def library_model(mosaic_model):
    model2 = mosaic_model.copy()
    model2.meta.instrument.optical_element = "F158"
    return ModelLibrary([mosaic_model, model2])


@pytest.fixture
def library_model_all_nan(mosaic_model):
    model1 = mosaic_model.copy()
    model1.data[:] = np.nan
    model2 = mosaic_model.copy()
    model2.data[:] = np.nan
    model2.meta.instrument.optical_element = "F158"
    return ModelLibrary([model1, model2])


def check_psf_matched_catalog(cat, all_filters, ref_filter):
    """
    Assert that a multiband catalog has correct PSF-matched column structure.

    Parameters
    ----------
    cat : `~astropy.table.Table`
    all_filters : list of str
        Lowercase filter names present in the catalog (e.g. ["f158", "f184"]).
    ref_filter : str
        Lowercase reference filter (e.g. "f184").
    """
    n_aper = len(cat.meta["aperture_radii"]["circle_pix"])
    matched_bands = [f"{f}m" for f in all_filters if f != ref_filter]

    assert cat.meta.get("psf_match_reference_filter") == ref_filter.upper()

    # All original filter bands should have aperture columns
    for f in all_filters:
        assert (
            sum(1 for c in cat.colnames if match(rf"^aper\d+_{f}_flux$", c)) == n_aper
        )

    for mb in matched_bands:
        # Aperture and background columns must be present
        assert (
            sum(1 for c in cat.colnames if match(rf"^aper\d+_{mb}_flux$", c)) == n_aper
        )
        assert f"aper_bkg_{mb}_flux" in cat.colnames
        assert f"aper_bkg_{mb}_flux_err" in cat.colnames
        # kron and segment flux columns must be present
        assert f"kron_{mb}_flux" in cat.colnames
        assert f"kron_{mb}_flux_err" in cat.colnames
        assert f"segment_{mb}_flux" in cat.colnames
        assert f"segment_{mb}_flux_err" in cat.colnames
        # abmag columns must not be present (redundant; derivable from flux)
        assert not any("_abmag" in c and mb in c for c in cat.colnames)
        # PSF flux columns must not be present (PSFs are already matched to each filter)
        assert f"psf_{mb}_flux" not in cat.colnames
        assert f"psf_{mb}_flux_err" not in cat.colnames
        # othershape columns (sharpness etc.) not computed for matched bands
        for param in ["sharpness", "roundness1", "is_extended", "fluxfrac_radius_50"]:
            assert f"{param}_{mb}" not in cat.colnames

    # Reference filter must have no matched aperture columns
    assert (
        sum(1 for c in cat.colnames if match(rf"^aper\d+_{ref_filter}m_flux$", c)) == 0
    )

    # All matched bands must have the same column structure
    if len(matched_bands) > 1:
        template = sorted(
            c.replace(matched_bands[0], "BAND")
            for c in cat.colnames
            if matched_bands[0] in c
        )
        for mb in matched_bands[1:]:
            assert (
                sorted(c.replace(mb, "BAND") for c in cat.colnames if mb in c)
                == template
            )

    # ee_fractions: exactly the original filter keys, correct lengths
    assert "ee_fractions" in cat.meta
    assert isinstance(cat.meta["ee_fractions"], dict)
    assert set(cat.meta["ee_fractions"].keys()) == set(all_filters)
    for value in cat.meta["ee_fractions"].values():
        assert len(value) == n_aper


def shared_tests(
    result, cat, library_model, save_results, function_jail, shape=(101, 101)
):
    with library_model:
        input_model = library_model.borrow(0)
        assert result.meta.data_release_id == input_model.meta.data_release_id
        library_model.shelve(input_model, modify=False)

    check_psf_matched_catalog(cat, ["f158", "f184"], "f184")

    if len(cat) > 0:
        assert np.min(cat["x_centroid"]) > 0.0
        assert np.min(cat["y_centroid"]) > 0.0
        assert np.max(cat["x_centroid"]) < shape[0]
        assert np.max(cat["y_centroid"]) < shape[1]

        for colname in cat.colnames:
            if (
                "flux" in colname
                and "fluxfrac" not in colname
                and "aper_bkg" not in colname
            ):
                assert cat[colname].unit == "nJy"
                assert "f158" in colname or "f184" in colname
            if colname.endswith("_flux"):
                assert f"{colname}_err" in cat.colnames

    catalog_filepath = Path(function_jail / f"{result.meta.filename}_cat.parquet")
    segmentation_map_filepath = Path(
        function_jail / f"{result.meta.filename}_segm.asdf"
    )

    if save_results:
        assert catalog_filepath.exists()
        tbl = pyarrow.parquet.read_table(catalog_filepath)
        assert isinstance(tbl, pyarrow.Table)

        assert segmentation_map_filepath.exists()
        segm_model = rdm.open(segmentation_map_filepath)
        assert isinstance(segm_model, MultibandSegmentationMapModel)
        assert (
            segm_model.meta.get("psf_match_reference_filter")
            == cat.meta["psf_match_reference_filter"]
        )
    else:
        assert not catalog_filepath.exists()
        assert not segmentation_map_filepath.exists()


@pytest.mark.parametrize("fit_psf", (True, False))
@pytest.mark.parametrize(
    "snr_threshold, npixels, save_results",
    (
        (3, 10, True),
        (7, 10, False),
    ),
)
def test_multiband_catalog(
    library_model, fit_psf, snr_threshold, npixels, save_results, function_jail
):
    result, _ = MultibandCatalogStep.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=snr_threshold,
        npixels=npixels,
        fit_psf=fit_psf,
        save_results=save_results,
        deblend=True,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7

    shared_tests(result, cat, library_model, save_results, function_jail)


def test_multiband_catalog_populates_dust_ebv(library_model, function_jail):
    """Ensure the joined multiband catalog contains the detection-level dust_ebv."""
    result, _ = MultibandCatalogStep.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=False,
        save_results=False,
        deblend=True,
    )
    cat = result.source_catalog
    assert "dust_ebv" in cat.colnames
    assert len(cat["dust_ebv"]) == len(cat)
    assert cat["dust_ebv"].dtype == np.float32


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_no_detections(library_model, save_results, function_jail):
    result, _ = MultibandCatalogStep.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=1000,  # high threshold to ensure no detections
        npixels=10,
        fit_psf=False,
        save_results=save_results,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 0


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_invalid_inputs(
    library_model_all_nan, save_results, function_jail
):
    result, _ = MultibandCatalogStep.call(
        library_model_all_nan,
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=False,
        save_results=save_results,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 0


@pytest.mark.parametrize("save_results", (True, False))
def test_multiband_catalog_some_invalid_inputs(
    library_model, save_results, function_jail
):
    # Modify the first model in the library to have all NaN values
    with library_model:
        model = library_model.borrow(0)  # f184 model
        model.data[:] = np.nan
        model.err[:] = np.nan
        model.var_rnoise[:] = np.nan
        library_model.shelve(model, modify=True)

    result, _ = MultibandCatalogStep.call(
        library_model,
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=False,
        save_results=save_results,
        deblend=True,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7
    assert "segment_f158_flux" in cat.colnames
    assert "segment_f184_flux" in cat.colnames
    assert np.all(np.isnan(cat["segment_f184_flux"]))
    assert np.all(np.isnan(cat["segment_f184_flux_err"]))


def mosaic_si_nan_model(shape=(500, 500)):
    model = MosaicModel.create_fake_data(shape=shape)
    data, err = make_si_test_image()
    model.data = data
    model.err = err
    model.var_rnoise = err**2
    model.var_poisson = err**2
    model.weight = 1.0 / err
    model.meta.instrument.optical_element = "F184"
    model.meta.coadd_info.time_first = Time("2027-01-01T00:00:00")
    model.meta.wcsinfo.pixel_scale = 0.11 / 3600  # degrees

    model.meta.wcsinfo.ra_ref = 270.0  # degrees
    model.meta.wcsinfo.dec_ref = 66.0  # degrees
    model.meta.wcsinfo.roll_ref = 0.0  # degrees
    model.meta.coadd_info.exposure_time = EXPTIME  # seconds

    model.meta.resample.pixfrac = 0.5
    model.meta.data_release_id = "r1"

    # Create WCS
    model.meta.wcs = mk_gwcs(
        model.meta.wcsinfo.ra_ref,
        model.meta.wcsinfo.dec_ref,
        model.meta.wcsinfo.roll_ref,
        bounding_box=((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5)),
        shape=shape,
    )

    return model


@pytest.fixture
def library_model2():
    si_model1 = mosaic_si_nan_model()
    si_model2 = deepcopy(si_model1)
    si_model2.meta.instrument.optical_element = "F158"
    return ModelLibrary([si_model1, si_model2])


@pytest.fixture
def libraries_si_nan():
    libs = {}
    nan_y_pos = nan_x_pos = None

    model_base = mosaic_si_nan_model()

    for si_type in ['NoNan', 'Grid', 'Block']:
        model1 = deepcopy(model_base)
        model1.meta.instrument.optical_element = "F158"

        if nan_y_pos is None:
            nan_y_pos, nan_x_pos = make_source_grid(
                model1,
                yxmax=model1.data.shape,
                yxoffset=(50, 50),
                yxgrid=(20, 20),
                seed=RNG_SEED,
            )
            nan_y_pos, nan_x_pos = np.round(nan_y_pos).astype(int), np.round(nan_x_pos).astype(int)

        # NaN in a quadrant
        if si_type == 'Block':
            model1.data[250:, 250:] = np.nan

        # NaN on grid diagonal
        if si_type != 'NoNan':
            model1.data[nan_y_pos[0::21], nan_x_pos[0::21]] = np.nan
            nanmask = np.isnan(model1.data)

        nanmask = np.isnan(model1.data)

        # Second model
        model2 = deepcopy(model1)
        model2.meta.instrument.optical_element = "F184"

        nanmask = np.isnan(model2.data)

        libs[si_type] = ModelLibrary([model1, model2])

    return libs


@pytest.mark.parametrize("fit_psf", (True, False))
@pytest.mark.parametrize(
    "snr_threshold, npixels, save_results",
    (
        (7, 5, False),
        (15, 4, True),
    ),
)
def test_multiband_source_injection_catalog(
    library_model2, fit_psf, snr_threshold, npixels, save_results, function_jail
):
    result, _ = MultibandCatalogStep.call(
        library_model2,
        bkg_boxsize=30,
        snr_threshold=snr_threshold,
        npixels=npixels,
        fit_psf=fit_psf,
        deblend=True,
        inject_sources=True,
        inject_seed=50,
        save_results=save_results,
        save_debug_info=True,
    )

    # Original objects
    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 112

    # Source injected and original images
    si_cat = result.source_injection_catalog
    assert isinstance(si_cat, Table)
    assert len(si_cat) >= len(cat)

    if save_results:
        filepath = Path(function_jail / f"{result.meta.filename}_segm.asdf")
        assert filepath.exists()
        segm_mod = rdm.open(filepath)
        assert isinstance(segm_mod, MultibandSegmentationMapModel)
        assert segm_mod.data.shape == segm_mod.si_data.shape
        assert isinstance(segm_mod.injected_sources, Table)
        assert len(segm_mod.injected_sources[0]) <= len(si_cat)

        assert np.count_nonzero(
            segm_mod.recovered_sources["best_injected_index"] != -1
        ) >= (400 / 2)

    # Old lines from other MBC tests
    shared_tests(
        result,
        cat,
        library_model2,
        save_results,
        function_jail,
        shape=(5000, 5000),
    )


def test_multiband_source_injection_nan_catalog(libraries_si_nan, function_jail):
    step = MultibandCatalogStep()

    # Dictionaries to hold step output
    res_cat = {}
    results = {}

    # Create three sets of libraries with two copies of the same mosaic:
    # NoNan - original mosaics
    # Grid - NoNan mosaics with the center locations of injected
    #     sources on a diagonal set to NaN
    # Block - Grid mosaics with a quadrant also set to NaN
    libmods = libraries_si_nan

    # Obtain source injection locations
    with libmods['NoNan']:
        for si_model in libmods['NoNan']:
            si_y_pos, si_x_pos = make_source_grid(
                        si_model,
                        yxmax=si_model.data.shape,
                        yxoffset=(50, 50),
                        yxgrid=(20, 20),
                        seed=RNG_SEED,
                    )
            y_pos_idx, x_pos_idx = np.round(si_y_pos).astype(int), np.round(si_x_pos).astype(int)

            libmods['NoNan'].shelve(si_model, modify=False)

    # Specify the Grid SI locations that should have NaN pixels
    # From above: yxgrid=(20, 20), so every 21st entry is a diagonal
    y_nan_idx, x_nan_idx = y_pos_idx[0::21], x_pos_idx[0::21]

    # Run the MultibandCatalogStep on all three libraries
    for si_type in ['NoNan', 'Grid', 'Block']:
        res_cat[si_type], results[si_type] = step.call(
            libmods[si_type],
            bkg_boxsize=30,
            snr_threshold=15,
            npixels=4,
            fit_psf=False,
            deblend=True,
            inject_sources=True,
            inject_seed=RNG_SEED,
            save_results=False,
            save_debug_info=True,
        )

    # Create masks for data comparisons across mosaics

    # Grid Nan mask
    grid_nan_stamp_mask = np.zeros_like(results["Grid"].si_detection_image, dtype=bool)
    for x,y in zip(y_nan_idx, x_nan_idx):
        grid_nan_stamp_mask[y-10:y+11, x-10:x+11] = True

    # Block Nan mask
    block_nan_stamp_mask = deepcopy(grid_nan_stamp_mask)
    block_nan_stamp_mask[240:, 240:] = True

    # Compare data

    # Grid vs NoNan
    # All pixels far from the NaN injected sources should be close
    assert np.allclose(results["NoNan"].si_detection_image[~grid_nan_stamp_mask],
       results["Grid"].si_detection_image[~grid_nan_stamp_mask], atol=MEANFLUX * EXPTIME)

    # Block vs NoNan
    # All pixels far from the NaN injected sources and outside of
    # the large NaN region should be close
    assert np.allclose(results["NoNan"].si_detection_image[~block_nan_stamp_mask],
       results["Block"].si_detection_image[~block_nan_stamp_mask], atol=MEANFLUX * EXPTIME)

    # Catalog tests to ensure that the sources injected at
    # NaN pixel center locations are reasonably close

    # Create catalogs of "best" selected recovered sources (only one source per injection)
    nonan_rs_no_dlbs = results['NoNan'].recovered_sources[results['NoNan'].recovered_sources['best_injected_index'] != -1]
    grid_rs_no_dlbs = results['Grid'].recovered_sources[results['Grid'].recovered_sources['best_injected_index'] != -1]
    block_rs_no_dlbs = results['Block'].recovered_sources[results['Block'].recovered_sources['best_injected_index'] != -1]

    # Create intermediate lists of NaN pixel indices on each mosaic group
    nan_idx = [i * 21 for i in range(0, 20)]
    nonan_ce = list(set(nan_idx) & set(nonan_rs_no_dlbs['best_injected_index'].tolist()))
    grid_ce = list(set(nan_idx) & set(grid_rs_no_dlbs['best_injected_index'].tolist()))
    block_ce = list(set(nan_idx) & set(block_rs_no_dlbs['best_injected_index'].tolist()))

    # Create matched lists of NaN pixel indices between
    # NoNan & Grid and NoNan & Block
    grid_inter_ce = list(set(nonan_ce) & set(grid_ce))
    block_inter_ce = list(set(nonan_ce) & set(block_ce))

    # Create numpy indices for each matched catalog object
    ngce_idx = np.array([np.where(nonan_rs_no_dlbs['best_injected_index'] == ce)[0][0] for ce in grid_inter_ce])
    gce_idx = np.array([np.where(grid_rs_no_dlbs['best_injected_index'] == ce)[0][0] for ce in grid_inter_ce])
    nbce_idx = np.array([np.where(nonan_rs_no_dlbs['best_injected_index'] == ce)[0][0] for ce in block_inter_ce])
    bce_idx = np.array([np.where(block_rs_no_dlbs['best_injected_index'] == ce)[0][0] for ce in block_inter_ce])

    # Remove obvious object mismatch
    # (SI can have multiple matches for one object, and it chooses the closest.
    #  NaNs can cause a mismatch.)
    bad_idx_mask = np.isclose(grid_rs_no_dlbs['ellipticity'][gce_idx],
                            nonan_rs_no_dlbs['ellipticity'][ngce_idx], atol=0.3)
    ngce_idx = ngce_idx[bad_idx_mask]
    gce_idx = gce_idx[bad_idx_mask]

    # Grid catalog tests

    # Ensure Grid has NaN flux within a 0.4" aperture of NaN pixels
    assert np.all(np.isnan(grid_rs_no_dlbs['aper04_f158_flux'][gce_idx]))

    # Ensure Grid Nan pixel objects have similar locations to matched NoNan object locations
    assert np.allclose(grid_rs_no_dlbs['x_centroid'][gce_idx],nonan_rs_no_dlbs['x_centroid'][ngce_idx], atol=2)
    assert np.allclose(grid_rs_no_dlbs['y_centroid'][gce_idx],nonan_rs_no_dlbs['y_centroid'][ngce_idx], atol=2)

    # Ensure Grid Nan pixel objects have similar size to matched NoNan object size
    assert np.allclose(grid_rs_no_dlbs['fluxfrac_radius_50_f184'][gce_idx],
                    nonan_rs_no_dlbs['fluxfrac_radius_50_f184'][ngce_idx],
                    rtol=0.4)

    # Ensure Grid Nan pixel objects have similar kron AB Mags to matched NoNan object kron AB Mags
    assert np.allclose(grid_rs_no_dlbs['kron_f184_abmag'][gce_idx],
                    nonan_rs_no_dlbs['kron_f184_abmag'][ngce_idx],
                    rtol=0.2)

    # Block catalog tests

    # Ensure Block has NaN flux within a 0.4" aperture of NaN pixels
    assert np.all(np.isnan(block_rs_no_dlbs['aper04_f158_flux'][bce_idx]))

    # Ensure Block Nan pixel objects have similar locations to matched NoNan object locations
    assert np.allclose(block_rs_no_dlbs['x_centroid'][bce_idx],nonan_rs_no_dlbs['x_centroid'][nbce_idx], atol=2)
    assert np.allclose(block_rs_no_dlbs['y_centroid'][bce_idx],nonan_rs_no_dlbs['y_centroid'][nbce_idx], atol=2)

    # Ensure Block Nan pixel objects have similar size to matched NoNan object size
    assert np.allclose(block_rs_no_dlbs['fluxfrac_radius_50_f184'][bce_idx],
                    nonan_rs_no_dlbs['fluxfrac_radius_50_f184'][nbce_idx],
                    rtol=0.4)

    # Ensure Block Nan pixel objects have similar kron AB Mags to matched NoNan object kron AB Mags
    assert np.allclose(block_rs_no_dlbs['kron_f184_abmag'][bce_idx],
                    nonan_rs_no_dlbs['kron_f184_abmag'][nbce_idx],
                    rtol=0.2)


def test_match_recovered_sources():
    # Create minimal source catalog data
    shape = (2000, 2000)
    orig_xy = [[510, 510], [400, 400], [1000, 1000], [450, 1450], [100, 100]]
    injected_xy = [[500, 500], [500, 1500], [1500, 500], [1500, 1500]]
    injected_hlr = [2.0, 12.0, 5.0, 1.0]
    si_cat_xy = [
        [510.1, 510.1],
        [400.1, 400.1],
        [1000.1, 1000.1],
        [450.1, 1450.1],
        [100.1, 100.1],
        [500.1, 500.1],
        [500.1, 1500.1],
        [1500.1, 500.1],
        [1500.1, 1500.1],
    ]

    # Make wcs object
    wcsobj = mk_gwcs(
        270.0,  # degrees
        66.0,  # degrees
        0.0,  # degrees
        bounding_box=((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5)),
        shape=shape,
    )

    # Create source catalogs
    orig_table = Table()
    injected_table = Table()
    si_cat_table = Table()

    # Conert x, y to ra & dec
    orig_table["x"], orig_table["y"] = zip(*orig_xy, strict=True)
    orig_table["ra"], orig_table["dec"] = wcsobj.pixel_to_world_values(
        np.array(orig_table["x"]), np.array(orig_table["y"])
    )

    injected_table["x"], injected_table["y"] = zip(*injected_xy, strict=True)
    injected_table["ra"], injected_table["dec"] = wcsobj.pixel_to_world_values(
        np.array(injected_table["x"]), np.array(injected_table["y"])
    )

    si_cat_table["x"], si_cat_table["y"] = zip(*si_cat_xy, strict=True)
    si_cat_table["ra"], si_cat_table["dec"] = wcsobj.pixel_to_world_values(
        np.array(si_cat_table["x"]), np.array(si_cat_table["y"])
    )

    # Additional table columns
    injected_table["half_light_radius"] = injected_hlr
    injected_table["half_light_radius"].unit = "arcsec"
    orig_table["empty"] = 0
    si_cat_table["one"] = 1

    # Conduct matching
    rec_table = match_recovered_sources(orig_table, injected_table, si_cat_table)

    # Test that distant sources were dropped
    assert len(rec_table) < len(si_cat_table)

    # Test that grid sources were found
    assert np.count_nonzero(rec_table["best_injected_index"] != -1) == len(
        injected_table
    )

    # Test columns included or excluded as expected
    assert "one" in rec_table.colnames
    assert "empty" not in rec_table.colnames


@pytest.fixture
def library_model_f062_f129_f213(mosaic_model):
    """
    Library with F062, F129, F213 for column content tests.
    """
    model1 = mosaic_model.copy()
    model1.meta.instrument.optical_element = "F062"

    model2 = mosaic_model.copy()
    model2.meta.instrument.optical_element = "F129"

    model3 = mosaic_model.copy()
    model3.meta.instrument.optical_element = "F213"

    return ModelLibrary([model1, model2, model3])


@pytest.fixture
def library_model_three_filters(mosaic_model):
    """
    Library with F062, F158, F184.
    """
    model1 = mosaic_model.copy()
    model1.meta.instrument.optical_element = "F062"

    model2 = mosaic_model.copy()
    model2.meta.instrument.optical_element = "F158"

    model3 = mosaic_model.copy()
    model3.meta.instrument.optical_element = "F184"

    # input models not in wavelength order to test sorting
    return ModelLibrary([model2, model3, model1])


@pytest.mark.parametrize("fit_psf", (True, False))
@pytest.mark.parametrize(
    "psf_match_reference_filter",
    [
        None,  # reddest (F184 auto-selected)
        "F158",  # middle
        "F062",  # bluest
    ],
)
def test_multiband_catalog_reference_filter(
    library_model_three_filters,
    fit_psf,
    psf_match_reference_filter,
    function_jail,
):
    """
    Test PSF matching with reddest, middle, and bluest reference filters.

    Uses [F062, F158, F184]; parametrized over which filter is the reference.
    Non-reference filters should have PSF-matched columns; the reference
    should not.  ee_fractions should contain the three original filter keys
    only (no 'm' variants).
    """
    kwargs = dict(
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=fit_psf,
        deblend=True,
        save_results=False,
    )
    if psf_match_reference_filter is not None:
        kwargs["psf_match_reference_filter"] = psf_match_reference_filter

    result, _ = MultibandCatalogStep.call(library_model_three_filters, **kwargs)

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7

    # None means reddest (F184) is auto-selected
    ref = (psf_match_reference_filter or "F184").lower()
    check_psf_matched_catalog(cat, ["f062", "f158", "f184"], ref)


@pytest.mark.parametrize("fit_psf", (True, False))
def test_multiband_catalog_column_content(
    library_model_f062_f129_f213, fit_psf, function_jail
):
    """
    Test that matched catalog columns contain the fields we need.

    With F062, F129, F213 and F129 as the reference:
    - F062m: bluer than reference, normal PSF convolution
    - F213m: redder than reference, synthetic correction factors
    - F129 (reference): only original measurements, no matched columns.
    """
    result, _ = MultibandCatalogStep.call(
        library_model_f062_f129_f213,
        bkg_boxsize=50,
        snr_threshold=3,
        npixels=10,
        fit_psf=fit_psf,
        deblend=True,
        psf_match_reference_filter="F129",
        save_results=False,
    )

    cat = result.source_catalog
    assert isinstance(cat, Table)
    assert len(cat) == 7

    check_psf_matched_catalog(cat, ["f062", "f129", "f213"], "f129")
