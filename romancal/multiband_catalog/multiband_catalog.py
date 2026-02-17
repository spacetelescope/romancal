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
from romancal.multiband_catalog.detection_image import make_detection_image
from romancal.multiband_catalog.utils import add_filter_to_colnames
from romancal.source_catalog import injection
from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import make_segmentation_image
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.utils import get_ee_spline

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


_SKIP_IMAGE_META_KEYS = {"wcs", "individual_image_meta"}
_SKIP_BLEND_KEYS = {"wcsinfo"}


def multiband_catalog(self, library, example_model, cat_model, ee_spline):
    """
    Create a multiband catalog of sources including photometry and basic
    shape measurements.

    Parameters
    -----------
    library : `~romancal.datamodels.ModelLibrary`
        The library of models.
    example_model : `MosaicModel` or `ImageModel`
        Example model.
    cat_model : `MultibandSourceCatalogModel`
        Catalog model.
    ee_spline : `astropy.modeling.fitting.SplineSplrepFitter

    Returns
    -------
    segment_img : `SegmentationImage` or set
        The segmentation image.
    cat_model : `MultibandSourceCatalogModel`
        Updated catalog.
    msg : str (optional)
        Reason for empty file.
    """
    # All input MosaicImages in the ModelLibrary are assumed to have
    # the same shape and be pixel aligned.

    log.info("Calculating and subtracting background")
    library = subtract_background_library(library, self.bkg_boxsize)

    log.info("Creating detection image")
    # Define the kernel FWHMs for the detection image
    # TODO: sensible defaults
    # TODO: redefine in terms of intrinsic FWHM
    if self.kernel_fwhms is None:
        self.kernel_fwhms = [2.0, 5.0]

    # TODO: det_img is saved in the MosaicSegmentationMapModel;
    # do we also want to save the det_err?
    det_img = make_detection_image(library, self.kernel_fwhms)

    # Estimate background rms from detection image to calculate a
    # threshold for source detection
    mask = ~np.isfinite(det_img)

    # Return an empty segmentation image and catalog table if all
    # pixels are masked in the detection image.
    if np.all(mask):
        msg = (
            "Cannot create source catalog. All "
            "pixels in the detection image are masked."
        )
        return det_img.shape, cat_model, msg

    bkg = RomanBackground(
        det_img,
        box_size=self.bkg_boxsize,
        coverage_mask=mask,
    )
    bkg_rms = bkg.background_rms

    log.info("Detecting sources")
    segment_img = make_segmentation_image(
        det_img,
        snr_threshold=self.snr_threshold,
        npixels=self.npixels,
        bkg_rms=bkg_rms,
        deblend=self.deblend,
        mask=mask,
    )

    if segment_img is None:  # no sources found
        msg = "Cannot create source catalog. No sources were detected."
        return det_img.shape, cat_model, msg

    segment_img.detection_image = det_img.copy()

    # Define the detection image model
    det_model = datamodels.MosaicModel.create_minimal()
    det_model.data = det_img
    det_model.err = np.ones_like(det_img)

    # TODO: this is a temporary solution to get model attributes
    # currently needed in RomanSourceCatalog
    det_model.weight = example_model.weight
    det_model.meta = example_model.meta

    # The stellar FWHM is needed to define the kernel used for
    # the DAOStarFinder sharpness and roundness properties.
    # TODO: measure on a secondary detection image with minimal
    # smoothing?; use the same detection image for basic shape
    # measurements?
    star_kernel_fwhm = np.min(self.kernel_fwhms)

    log.info("Creating catalog for detection image")
    det_catobj = RomanSourceCatalog(
        det_model,
        cat_model,
        segment_img,
        det_img,
        star_kernel_fwhm,
        fit_psf=False,  # not needed for detection image
        detection_cat=None,
        mask=mask,
        cat_type="dr_det",
        ee_spline=ee_spline,
    )

    # Generate the catalog for the detection image.
    # We need to make this catalog before we pass det_catobj
    # to the RomanSourceCatalog constructor.
    det_cat = det_catobj.catalog
    det_cat.meta["ee_fractions"] = {}

    time_means = []
    exposure_times = []

    # Create catalogs for each input image
    with library:
        for model in library:
            mask = ~np.isfinite(model.data) | ~np.isfinite(model.err) | (model.err <= 0)

            if self.fit_psf:
                filter_name = model.meta.instrument.optical_element
                log.info(f"Creating catalog for {filter_name} image")
                ref_file = self.get_reference_file(model, "epsf")
                log.info("Using ePSF reference file: %s", ref_file)
                psf_ref_model = datamodels.open(ref_file)
            else:
                psf_ref_model = None

            apcorr_ref = self.get_reference_file(model, "apcorr")
            ee_spline = get_ee_spline(model, apcorr_ref)

            catobj = RomanSourceCatalog(
                model,
                cat_model,
                segment_img,
                None,
                star_kernel_fwhm,
                fit_psf=self.fit_psf,
                detection_cat=det_catobj,
                mask=mask,
                psf_ref_model=psf_ref_model,
                cat_type="dr_band",
                ee_spline=ee_spline,
            )

            # Add the filter name to the column names
            filter_name = model.meta.instrument.optical_element
            cat = add_filter_to_colnames(catobj.catalog, filter_name)
            ee_fractions = cat.meta["ee_fractions"]

            # TODO: what metadata do we want to keep, if any,
            # from the filter catalogs?
            cat.meta = None

            # Add the filter catalog to the multiband catalog.
            # The outer join prevents an empty table if any
            # columns have the same name but different values
            # (e.g., repeated filter names)
            det_cat = join(det_cat, cat, keys="label", join_type="outer")
            det_cat.meta["ee_fractions"][filter_name.lower()] = ee_fractions

            # accumulate image metadata
            image_meta = {
                k: copy.deepcopy(v)
                for k, v in model["meta"].items()
                if k not in _SKIP_IMAGE_META_KEYS
            }
            cat_model.meta.image_metas.append(image_meta)

            # blend model with catalog metadata
            if model.meta.file_date < cat_model.meta.image.file_date:
                cat_model.meta.image.file_date = model.meta.file_date

            for key, value in image_meta.items():
                if key in _SKIP_BLEND_KEYS:
                    continue
                if not isinstance(value, dict):
                    # skip blending of single top-level values
                    continue
                if key not in cat_model.meta:
                    # skip blending if the key is not in the catalog meta
                    continue
                if key == "coadd_info":
                    cat_model.meta[key]["time_first"] = min(
                        cat_model.meta[key]["time_first"], value["time_first"]
                    )
                    cat_model.meta[key]["time_last"] = max(
                        cat_model.meta[key]["time_last"], value["time_last"]
                    )
                    time_means.append(value["time_mean"])
                    exposure_times.append(value["exposure_time"])
                else:
                    # set non-matching metadata values to None
                    for subkey, subvalue in value.items():
                        if cat_model.meta[key].get(subkey, None) != subvalue:
                            cat_model.meta[key][subkey] = None

            library.shelve(model, modify=False)

    # finish blending
    cat_model.meta.coadd_info.time_mean = Time(time_means).mean()
    cat_model.meta.coadd_info.exposure_time = np.mean(exposure_times)

    # Put the resulting multiband catalog in the model
    cat_model.source_catalog = det_cat

    return segment_img, cat_model, None


def make_source_injected_library(library, seed=None):
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
                # si_model.var_poisson = si_model.err**2
                # XXX UNDO
                # si_model.var_poisson = si_model.err**2
                si_model.var_poisson = np.zeros_like(si_model.err.shape)

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
            si_model = injection.inject_sources(si_model, si_cat, seed,)

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
