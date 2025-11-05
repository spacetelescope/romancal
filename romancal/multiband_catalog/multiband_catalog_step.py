"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING

import numpy as np
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
from romancal.source_catalog.save_utils import save_all_results, save_empty_results
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.utils import get_ee_spline
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["MultibandCatalogStep"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


_SKIP_IMAGE_META_KEYS = {"wcs", "individual_image_meta"}
_SKIP_BLEND_KEYS = {"wcsinfo"}


class MultibandCatalogStep(RomanStep):
    """
    Create a multiband catalog of sources including photometry and basic
    shape measurements.

    Parameters
    -----------
    input : str or `~romancal.datamodels.ModelLibrary`
        Path to an ASDF file or a `~romancal.datamodels.ModelLibrary`
        that contains `~roman_datamodels.datamodels.MosaicImageModel`
        models.
    """

    class_alias = "multiband_catalog"
    reference_file_types: ClassVar = ["apcorr"]

    spec = """
        bkg_boxsize = integer(default=100)   # background mesh box size in pixels
        kernel_fwhms = float_list(default=None)  # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # per-pixel SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        suffix = string(default='cat')        # Default suffix for output files
        fit_psf = boolean(default=True)       # fit source PSFs for accurate astrometry?
        inject_sources = boolean(default=False) # Inject sources into images
    """

    def process(self, library):
        # All input MosaicImages in the ModelLibrary are assumed to have
        # the same shape and be pixel aligned.
        if isinstance(library, str):
            library = ModelLibrary(library)
        if not isinstance(library, ModelLibrary):
            raise TypeError("library input must be a ModelLibrary object")

        with library:
            example_model = library.borrow(0)
            library.shelve(example_model, modify=False)

        # Initialize the source catalog model, copying the metadata
        # from the example model. Some of this may be overwritten
        # during metadata blending.
        cat_model = datamodels.MultibandSourceCatalogModel.create_minimal(
            {"meta": example_model.meta}
        )
        cat_model.meta["image"] = {
            # try to record association name else fall back to example model filename
            "filename": library.asn.get("table_name", example_model.meta.filename),
            "file_date": example_model.meta.file_date,
            # this may be overwritten during metadata blending
        }
        cat_model.meta["image_metas"] = []
        # copy over data_release_id, ideally this will come from the association
        if "data_release_id" in example_model.meta:
            cat_model.meta.data_release_id = example_model.meta.data_release_id

        log.info("Creating ee_fractions model for first image")
        apcorr_ref = self.get_reference_file(example_model, "apcorr")
        ee_spline = get_ee_spline(example_model, apcorr_ref)

        # Define the output filename for the source catalog model
        try:
            cat_model.meta.filename = library.asn["products"][0]["name"]
        except (AttributeError, KeyError):
            cat_model.meta.filename = "multiband_catalog"

        log.info("Calculating and subtracting background")
        library = subtract_background_library(library, self.bkg_boxsize)

        log.info("Creating detection image")
        # Define the kernel FWHMs for the detection image
        # TODO: sensible defaults
        # TODO: redefine in terms of intrinsic FWHM
        if self.kernel_fwhms is None:
            self.kernel_fwhms = [2.0, 20.0]

        # TODO: det_img is saved in the MosaicSegmentationMapModel;
        # do we also want to save the det_err?
        det_img, det_err = make_detection_image(library, self.kernel_fwhms)

        # Estimate background rms from detection image to calculate a
        # threshold for source detection
        mask = ~np.isfinite(det_img) | ~np.isfinite(det_err) | (det_err <= 0)

        # Return an empty segmentation image and catalog table if all
        # pixels are masked in the detection image.
        if np.all(mask):
            msg = "Cannot create source catalog. All pixels in the detection image are masked."
            return save_empty_results(self, det_img.shape, cat_model, msg=msg)

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
            return save_empty_results(self, det_img.shape, cat_model, msg=msg)

        segment_img.detection_image = det_img

        # Define the detection image model
        det_model = datamodels.MosaicModel.create_minimal()
        det_model.data = det_img
        det_model.err = det_err

        # TODO: this is a temporary solution to get model attributes
        # currently needed in RomanSourceCatalog
        det_model.weight = example_model.weight
        det_model.meta = example_model.meta

        # Source Injection
        if self.inject_sources:
            # Imports
            from astropy import units as u
            from astropy.coordinates import SkyCoord

            from romancal.skycell.tests.test_skycell_match import mk_gwcs

            # Make copy of detection image to inject sources into
            si_model = copy.deepcopy(det_model)

        # The stellar FWHM is needed to define the kernel used for
        # the DAOStarFinder sharpness and roundness properties.
        # TODO: measure on a secondary detection image with minimal
        # smoothing?; use the same detection image for basic shape
        # measurements?
        star_kernel_fwhm = np.min(self.kernel_fwhms)

        log.info("Creating catalog for detection image")
        det_catobj = RomanSourceCatalog(
            det_model,
            segment_img,
            det_img,
            star_kernel_fwhm,
            fit_psf=self.fit_psf,
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

        # Source Injection
        if self.inject_sources:
            # Generate SI background rms
            si_bkg = RomanBackground(
                si_model.data,
                box_size=self.bkg_boxsize,
                coverage_mask=mask,
            )
            si_bkg_rms = si_bkg.background_rms

            # Poisson variance required for source injection
            if "var_poisson" not in si_model:
                si_model.var_poisson = si_model.err**2

            # Create source grid points
            si_x_pos, si_y_pos = injection.make_source_grid(
                si_model, yxmax=si_model.data.shape, yxoffset=(50, 50), yxgrid=(20, 20)
            )

            si_cen = SkyCoord(
                ra=det_model.meta.wcsinfo.ra_ref * u.deg,
                dec=det_model.meta.wcsinfo.dec_ref * u.deg,
            )

            # Obtain exposure times and filters
            # This code assumes all filters have been coadded already,
            # and thus there is one image per filter
            si_filters = []
            si_exptimes = {}
            with library:
                for model in library:
                    si_filter_name = model.meta.instrument.optical_element

                    si_exptimes[si_filter_name] = float(
                        model.meta.coadd_info.exposure_time
                    )
                    si_filters.append(si_filter_name)

                    library.shelve(model, modify=False)

            # WCS object for ra & dec conversion
            wcsobj = mk_gwcs(
                det_model.meta.wcsinfo.ra_ref,
                det_model.meta.wcsinfo.dec_ref,
                det_model.meta.wcsinfo.roll_ref,
                bounding_box=(
                    (-0.5, model.data.shape[0] - 0.5),
                    (-0.5, model.data.shape[1] - 0.5),
                ),
                shape=model.data.shape,
            )

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

            # Save catalog to segmentation image
            segment_img.injected_sources = si_cat.as_array()

            # Inject sources into the detection image
            si_model = injection.inject_sources(si_model, si_cat)

            # Create SI segmentation image
            si_segment_img = make_segmentation_image(
                si_model.data,
                snr_threshold=self.snr_threshold,
                npixels=self.npixels,
                bkg_rms=si_bkg_rms,
                deblend=self.deblend,
                mask=mask,
            )

            # Save SI segmentation image
            si_segment_img.detection_image = si_model.data
            segment_img.si_detection_image = si_model.data

            log.info("Creating catalog for source injected detection image")
            si_det_catobj = RomanSourceCatalog(
                si_model,
                si_segment_img,
                si_model.data,
                star_kernel_fwhm,
                fit_psf=self.fit_psf,
                detection_cat=None,
                mask=mask,
                cat_type="dr_det",
                ee_spline=ee_spline,
            )

            si_det_cat = si_det_catobj.catalog
            si_det_cat.meta["ee_fractions"] = {}

            segment_img.si_segment_img = si_segment_img

        time_means = []
        exposure_times = []

        # Create catalogs for each input image
        with library:
            for model in library:
                if self.inject_sources:
                    si_model = copy.deepcopy(model)

                mask = (
                    ~np.isfinite(model.data)
                    | ~np.isfinite(model.err)
                    | (model.err <= 0)
                )

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

                if self.inject_sources:
                    # Poisson variance required for source injection
                    if "var_poisson" not in si_model:
                        si_model.var_poisson = si_model.err**2

                    # Inject sources into the model image
                    si_model = injection.inject_sources(si_model, si_cat)

                    # Create catalog
                    si_catobj = RomanSourceCatalog(
                        si_model,
                        si_segment_img,
                        None,
                        star_kernel_fwhm,
                        fit_psf=self.fit_psf,
                        detection_cat=si_det_catobj,
                        mask=mask,
                        psf_ref_model=psf_ref_model,
                        cat_type="dr_band",
                        ee_spline=ee_spline,
                    )

                    # Add the filter name to the column names
                    si_filter_name = si_model.meta.instrument.optical_element
                    si_model_cat = add_filter_to_colnames(
                        si_catobj.catalog, si_filter_name
                    )
                    si_ee_fractions = si_model_cat.meta["ee_fractions"]
                    si_model_cat.meta = None

                    # Add the filter catalog to the multiband catalog.
                    # The outer join prevents an empty table if any
                    # columns have the same name but different values
                    # (e.g., repeated filter names)
                    si_det_cat = join(
                        si_det_cat, si_model_cat, keys="label", join_type="outer"
                    )
                    si_det_cat.meta["ee_fractions"][si_filter_name.lower()] = (
                        si_ee_fractions
                    )

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

        if self.inject_sources:
            # Put the source injected multiband catalog in the model
            cat_model.source_injection_catalog = si_det_cat

        return save_all_results(self, segment_img, cat_model)
