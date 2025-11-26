"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import copy
import logging

import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.table import join
from astropy.time import Time
from roman_datamodels import datamodels

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog.background import subtract_background_library
from romancal.multiband_catalog.catalog_generator import create_filter_catalog
from romancal.multiband_catalog.detection_image import make_detection_image
from romancal.multiband_catalog.metadata import blend_image_metadata
from romancal.source_catalog import injection
from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import make_segmentation_image
from romancal.source_catalog.psf_matching import get_filter_wavelength
from romancal.source_catalog.source_catalog import RomanSourceCatalog

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def process_detection_image(self, library, example_model, ee_spline, catalog_model):
    """
    Create and process the detection image.

    This includes background estimation, source detection via
    segmentation, and creation of the detection catalog.

    Parameters
    ----------
    library : `romancal.datamodels.ModelLibrary`
        The library of models to process.

    example_model : `romancal.datamodels.MosaicImageModel`
        An example model from the library for metadata access.

    ee_spline : callable
        The encircled energy spline function.

    catalog_model : `romancal.datamodels.MultibandSourceCatalogModel`
        The output catalog model (for saving empty results if
        needed).

    Returns
    -------
    result : dict or tuple
        If successful, returns a dictionary with keys:
        - 'detection_model': The detection image model
        - 'mask': The total mask array
        - 'segment_img': The segmentation image
        - 'detection_catobj': The detection RomanSourceCatalog object
        - 'detection_catalog': The detection catalog table
        - 'star_kernel_fwhm': The stellar kernel FWHM

        If detection fails, returns the result of
        save_empty_results().
    """
    log.info("Creating detection image")

    # Define the kernel FWHMs for the detection image
    # TODO: sensible defaults
    # TODO: redefine in terms of intrinsic FWHM
    if self.kernel_fwhms is None:
        self.kernel_fwhms = [2.0, 5.0]

    detection_image = make_detection_image(library, self.kernel_fwhms)

    # Estimate background rms from detection image to calculate a
    # threshold for source detection
    mask = ~np.isfinite(detection_image)

    # Return image shape and empty catalog table if all pixels are
    # masked in the detection image
    if np.all(mask):
        msg = (
            "Cannot create source catalog. All "
            "pixels in the detection image are masked."
        )
        return detection_image.shape, catalog_model, msg

    log.info("Calculating background RMS for detection image")
    bkg = RomanBackground(
        detection_image,
        box_size=self.bkg_boxsize,
        coverage_mask=mask,
    )
    bkg_rms = bkg.background_rms

    log.info("Detecting sources")
    segment_img = make_segmentation_image(
        detection_image,
        snr_threshold=self.snr_threshold,
        npixels=self.npixels,
        bkg_rms=bkg_rms,
        deblend=self.deblend,
        mask=mask,
    )

    # Return image shape and empty catalog table if no sources
    # were detected
    if segment_img is None:  # no sources found
        msg = "Cannot create source catalog. No sources were detected."
        return detection_image.shape, catalog_model, msg

    segment_img.detection_image = detection_image.copy()

    # Define the detection image model
    detection_model = datamodels.MosaicModel.create_minimal()
    detection_model.data = detection_image
    detection_model.err = np.ones_like(detection_image)

    # TODO: this is a temporary solution to get model attributes
    # currently needed in RomanSourceCatalog
    detection_model.weight = example_model.weight
    detection_model.meta = example_model.meta

    # The stellar FWHM is needed to define the kernel used for
    # the DAOStarFinder sharpness and roundness properties.
    # TODO: measure on a secondary detection image with minimal
    # smoothing?; use the same detection image for basic shape
    # measurements?
    star_kernel_fwhm = np.min(self.kernel_fwhms)

    log.info("Creating catalog for detection image")
    detection_catobj = RomanSourceCatalog(
        detection_model,
        segment_img,
        detection_image,
        star_kernel_fwhm,
        fit_psf=False,  # not needed for detection image
        detection_cat=None,
        mask=mask,
        cat_type="dr_det",
        ee_spline=ee_spline,
    )

    # Generate the catalog for the detection image. The catalog
    # is lazily evalated, so we need to access it before we pass
    # detection_catobj to the RomanSourceCatalog constructor.
    detection_catalog = detection_catobj.catalog

    return {
        "detection_model": detection_model,
        "mask": mask,
        "segment_img": segment_img,
        "detection_catobj": detection_catobj,
        "detection_catalog": detection_catalog,
        "star_kernel_fwhm": star_kernel_fwhm,
    }


def finalize_ee_fractions(detection_catalog, filter_ee_fractions):
    """
    Consolidate and finalize encircled energy fractions.

    Accumulates ee_fractions from all filter processing and sorts
    them by wavelength in the detection catalog metadata. Only
    includes ee_fractions for original (non-PSF-matched) filter
    bands.

    The method modifies detection_catalog.meta["ee_fractions"] in
    place.

    Parameters
    ----------
    detection_catalog : Table
        The detection catalog where ee_fractions metadata will be
        stored.

    filter_ee_fractions : list of dict
        List of ee_fractions dictionaries from each filter,
        where each dict maps filter names to ee_fractions arrays.
    """
    detection_catalog.meta["ee_fractions"] = {}

    # Accumulate all ee_fractions from filters
    for ee_fracs in filter_ee_fractions:
        for key, value in ee_fracs.items():
            detection_catalog.meta["ee_fractions"][key] = value

    # Sort ee_fractions dictionary by filter wavelength
    if detection_catalog.meta.get("ee_fractions"):
        sorted_ee_fractions = dict(
            sorted(
                detection_catalog.meta["ee_fractions"].items(),
                key=lambda item: get_filter_wavelength(item[0]),
            )
        )
        detection_catalog.meta["ee_fractions"] = sorted_ee_fractions


def multiband_catalog(self, library, example_model, catalog_model, ee_spline):
    """
    Create a multiband catalog of sources including photometry and basic
    shape measurements.

    Parameters
    -----------
    library : `~romancal.datamodels.ModelLibrary`
        The library of models.
    example_model : `MosaicModel` or `ImageModel`
        Example model.
    catalog_model : `MultibandSourceCatalogModel`
        Catalog model.
    ee_spline : `astropy.modeling.fitting.SplineSplrepFitter

    Returns
    -------
    segment_img : `SegmentationImage` or set
        The segmentation image.
    catalog_model : `MultibandSourceCatalogModel`
        Updated catalog.
    msg : str (optional)
        Reason for empty file.

    Notes
    -----
    All input MosaicImages in the ModelLibrary are assumed to have the
    same shape and be pixel aligned.
    """
    log.info("Calculating and subtracting background")
    library = subtract_background_library(library, self.bkg_boxsize)

    # Create detection image, segmentation, and catalog
    detection_result = process_detection_image(
        self, library, example_model, ee_spline, catalog_model
    )

    # Check if detection failed (step returns save_empty_results)
    if isinstance(detection_result, tuple):
        log.warning("Detection image processing failed")
        return detection_result

    # Extract detection results
    segment_img = detection_result["segment_img"]
    detection_catobj = detection_result["detection_catobj"]
    detection_catalog = detection_result["detection_catalog"]
    star_kernel_fwhm = detection_result["star_kernel_fwhm"]

    time_means = []
    exposure_times = []
    filter_ee_fractions = []

    # Create catalogs for each input image
    with library:
        for model in library:
            # mask = ~np.isfinite(model.data) | ~np.isfinite(model.err) | (model.err <= 0)

            filter_name = model.meta.instrument.optical_element
            res = create_filter_catalog(
                model,
                filter_name,
                filter_name,  # ref_filter (same as filter_name to disable matching)
                0,  # ref_wavelength (ignored)
                segment_img,
                star_kernel_fwhm,
                detection_catobj,
                None,  # ref_model
                None,  # ref_filter_catalog
                None,  # ref_psf_model
                self.fit_psf,
                self.get_reference_file,
            )
            cat = res["catalog"]

            # TODO: what metadata do we want to keep, if any,
            # from the filter catalogs?
            cat.meta = None

            # Add the filter catalog to the multiband catalog.
            # The outer join prevents an empty table if any
            # columns have the same name but different values
            # (e.g., repeated filter names)
            detection_catalog = join(
                detection_catalog, cat, keys="label", join_type="outer"
            )

            # Store ee_fractions for consolidation
            filter_ee_fractions.append(res["ee_fractions"])

            # Accumulate and blend image metadata
            blend_image_metadata(model, catalog_model, time_means, exposure_times)

            library.shelve(model, modify=False)

    # Finish blending
    catalog_model.meta.coadd_info.time_mean = Time(time_means).mean()
    catalog_model.meta.coadd_info.exposure_time = np.mean(exposure_times)

    # Consolidate and sort ee_fractions
    finalize_ee_fractions(detection_catalog, filter_ee_fractions)

    # Put the resulting multiband catalog in the model
    catalog_model.source_catalog = detection_catalog

    return segment_img, catalog_model, None


def initialize_catalog_model(library, example_model):
    """
    Initialize the multiband source catalog model.

    Creates the catalog model and sets up initial metadata from
    the example model and association information.

    Parameters
    ----------
    library : `romancal.datamodels.ModelLibrary`
        The library of models to process.

    example_model : `romancal.datamodels.MosaicImageModel`
        An example model from the library for metadata access.

    Returns
    -------
    result : `romancal.datamodels.MultibandSourceCatalogModel`
        The initialized catalog model with set metadata.
    """
    # Initialize the source catalog model, copying the metadata
    # from the example model. Some of this may be overwritten
    # during metadata blending.
    catalog_model = datamodels.MultibandSourceCatalogModel.create_minimal(
        {"meta": example_model.meta}
    )
    catalog_model.meta["image"] = {
        # try to record association name else fall back to example model
        # filename
        "filename": library.asn.get("table_name", example_model.meta.filename),
        # file_date may be overwritten during metadata blending
        "file_date": example_model.meta.file_date,
    }
    catalog_model.meta["image_metas"] = []

    # copy over data_release_id, ideally this will come from the association
    if "data_release_id" in example_model.meta:
        catalog_model.meta.data_release_id = example_model.meta.data_release_id

    # Define the output filename for the source catalog model
    try:
        catalog_model.meta.filename = library.asn["products"][0]["name"]
    except (AttributeError, KeyError):
        catalog_model.meta.filename = "multiband_catalog"

    return catalog_model


def make_source_injected_library(library):
    """
    Create a library of source injected models.

    Parameters
    -----------
    input : str or `~romancal.datamodels.ModelLibrary`
        Path to an ASDF file or a `~romancal.datamodels.ModelLibrary`
        that contains `~roman_datamodels.datamodels.MosaicImageModel`
        models.

    Returns
    -------
    result : `~romancal.datamodels.ModelLibrary`
        The library of source injected models.

    si_cat : `~astropy.table.Table`
        Catalog of injected sources.
    """
    si_model_lst = []
    si_filters = []
    si_exptimes = {}
    si_cen = None

    # Cycle through library images to make source injected versions
    with library:
        for model in library:
            si_model = copy.deepcopy(model)
            library.shelve(model, modify=False)

            si_filter_name = si_model.meta.instrument.optical_element
            si_exptimes[si_filter_name] = float(si_model.meta.coadd_info.exposure_time)
            si_filters.append(si_filter_name)

            # Poisson variance required for source injection
            if "var_poisson" not in si_model:
                si_model.var_poisson = si_model.err**2

            # Set parameters for source injection
            # This only needs to be done once per library
            if si_cen is None:
                # Create source grid points
                si_x_pos, si_y_pos = injection.make_source_grid(
                    si_model,
                    yxmax=si_model.data.shape,
                    yxoffset=(50, 50),
                    yxgrid=(20, 20),
                )

                si_cen = coordinates.SkyCoord(
                    ra=si_model.meta.wcsinfo.ra_ref * u.deg,
                    dec=si_model.meta.wcsinfo.dec_ref * u.deg,
                )

                # Convert to RA & Dec
                wcsobj = si_model.meta.wcs
                si_ra, si_dec = wcsobj.pixel_to_world_values(
                    np.array(si_x_pos), np.array(si_y_pos)
                )

                # Generate cosmos-like catalog
                si_cat = injection.make_cosmoslike_catalog(
                    cen=si_cen,
                    ra=si_ra,
                    dec=si_dec,
                    exptimes=si_exptimes,
                )

                # Additional useful parameters
                si_cat["x_pos"] = si_x_pos
                si_cat["y_pos"] = si_y_pos
                si_cat["label"] = np.arange(len(si_x_pos))

            # Inject sources into the detection image
            si_model = injection.inject_sources(si_model, si_cat)

            # Add model to list for new library
            si_model_lst.append(si_model)

    # Return library of source injected models and injection catalog
    return ModelLibrary(si_model_lst), si_cat


def match_recovered_sources(original, injected, si_catalogs):
    """
    Create recovered sources which matches sources between
    an original catalog and injected catalog.

    Parameters
    -----------
    original : `~astropy.table.Table`
        Catalog of original sources.

    injected : `~astropy.table.Table`
        Catalog of injected sources.

    si_catalogs : `~astropy.table.Table`
        Catalog of sources after injection.

    Returns
    -------
    recovered : `~astropy.table.Table`
        Catalog of recovered sources.
    """

    # Create recovered catalog
    recovered = copy.deepcopy(si_catalogs)
    recovered["best_injected_index"] = -1

    def _make_skycoord(cat):
        ra, dec = cat["ra"], cat["dec"]
        if not hasattr(cat["ra"], "unit") or ra.unit is None:
            ra = ra * u.deg
            dec = dec * u.deg
        return coordinates.SkyCoord(ra, dec)

    # Create skycoord objects
    ocoord = _make_skycoord(original)
    icoord = _make_skycoord(injected)
    rcoord = _make_skycoord(recovered)

    # Set maximum radius for matching.
    # d = 3*sqrt(half_right_radius^2 + 0.2^2)
    # In addition, set a maximum of 10"
    max_sep_injected = 3 * np.sqrt(
        injected["half_light_radius"].to(u.arcsec).value ** 2 + 0.2**2
    )
    max_sep_injected = np.clip(max_sep_injected, 0, np.inf)

    # Trim recovered catalog to only objects near injected sources
    mi, mr, dist, _ = coordinates.search_around_sky(
        icoord, rcoord, np.max(max_sep_injected) * u.arcsec
    )
    m = dist < max_sep_injected[mi] * u.arcsec
    keep = np.zeros(len(rcoord), dtype="bool")
    keep[mr[m]] = True
    recovered = recovered[keep]
    rcoord = rcoord[keep]

    idx, sep, _ = coordinates.match_coordinates_sky(icoord, rcoord)
    m = sep < max_sep_injected * u.arcsec

    # Set matched injected sources
    recovered["best_injected_index"] = -1
    recovered["best_injected_index"][idx[m]] = np.flatnonzero(m)

    # Set distances
    idx, sep, _ = coordinates.match_coordinates_sky(rcoord, ocoord)
    recovered["dist_nearest"] = sep.to(u.arcsec)

    return recovered
