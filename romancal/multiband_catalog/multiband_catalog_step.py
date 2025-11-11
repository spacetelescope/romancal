"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from astropy.table import join
from astropy.time import Time
from roman_datamodels import datamodels

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog.background import subtract_background_library
from romancal.multiband_catalog.catalog_generator import create_filter_catalog
from romancal.multiband_catalog.detection_image import make_detection_image
from romancal.multiband_catalog.metadata import blend_image_metadata
from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import make_segmentation_image
from romancal.source_catalog.psf_matching import (
    get_filter_wavelength,
    get_reddest_filter,
)
from romancal.source_catalog.save_utils import save_all_results, save_empty_results
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.source_catalog.utils import get_ee_spline
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["MultibandCatalogStep"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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
        reference_filter = string(default=None)  # reference filter for PSF matching
    """

    def _initialize_catalog_model(self, library, example_model):
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
            # try to record association name else fall back to
            # example model filename
            "filename": library.asn.get("table_name", example_model.meta.filename),
            # this may be overwritten during metadata blending
            "file_date": example_model.meta.file_date,
        }
        catalog_model.meta["image_metas"] = []

        # copy over data_release_id, ideally this will come from the
        # association
        if "data_release_id" in example_model.meta:
            catalog_model.meta.data_release_id = example_model.meta.data_release_id

        # Define the output filename for the source catalog model
        try:
            catalog_model.meta.filename = library.asn["products"][0]["name"]
        except (AttributeError, KeyError):
            catalog_model.meta.filename = "multiband_catalog"

        return catalog_model

    def _process_detection_image(
        self, library, example_model, ee_spline, catalog_model
    ):
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
            self.kernel_fwhms = [2.0, 20.0]

        # TODO: detection_img is saved in the MosaicSegmentationMapModel;
        # do we also want to save the detection_err?
        detection_img, detection_err = make_detection_image(library, self.kernel_fwhms)

        # Estimate background rms from detection image to calculate a
        # threshold for source detection
        mask = (
            ~np.isfinite(detection_img)
            | ~np.isfinite(detection_err)
            | (detection_err <= 0)
        )

        # Return an empty segmentation image and catalog table if all
        # pixels are masked in the detection image.
        if np.all(mask):
            msg = (
                "Cannot create source catalog. All pixels in the "
                "detection image are masked."
            )
            return save_empty_results(self, detection_img.shape, catalog_model, msg=msg)

        log.info("Calculating background RMS for detection image")
        bkg = RomanBackground(
            detection_img,
            box_size=self.bkg_boxsize,
            coverage_mask=mask,
        )
        bkg_rms = bkg.background_rms

        log.info("Detecting sources")
        segment_img = make_segmentation_image(
            detection_img,
            snr_threshold=self.snr_threshold,
            npixels=self.npixels,
            bkg_rms=bkg_rms,
            deblend=self.deblend,
            mask=mask,
        )

        if segment_img is None:  # no sources found
            msg = "Cannot create source catalog. No sources were detected."
            return save_empty_results(self, detection_img.shape, catalog_model, msg=msg)

        segment_img.detection_image = detection_img

        # Define the detection image model
        detection_model = datamodels.MosaicModel()
        detection_model.data = detection_img
        detection_model.err = detection_err

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
            detection_img,
            star_kernel_fwhm,
            fit_psf=self.fit_psf,
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

    def _setup_reference_filter(self, library):
        """
        Determine and load the reference filter for PSF matching.

        Parameters
        ----------
        library : `romancal.datamodels.ModelLibrary`
            The library of models to process.

        Returns
        -------
        result : dict
            Dictionary with keys:
            - 'ref_filter': The reference filter name (uppercase)
            - 'ref_wavelength': The reference filter approx wavelength (nm)
            - 'ref_model': The reference filter model
            - 'ref_psf_model': The reference PSF model
        """
        # Determine reference filter for PSF matching
        if self.reference_filter is None:
            # Default to reddest filter
            ref_filter = get_reddest_filter(library)
            log.info(f"Using reddest filter as reference: {ref_filter}")
        else:
            # User specified reference filter
            ref_filter = self.reference_filter.upper()
            log.info(f"Using user-specified reference filter: {ref_filter}")

        # Load reference PSF model
        ref_model = None
        with library:
            for model in library:
                if model.meta.instrument.optical_element == ref_filter:
                    ref_model = model
                library.shelve(model, modify=False)
                if ref_model is not None:
                    break

        if ref_model is None:
            msg = (
                f"Reference filter {ref_filter} not found in library. "
                "Cannot perform PSF matching."
            )
            raise ValueError(msg)

        ref_psf_file = self.get_reference_file(ref_model, "epsf")
        ref_psf_model = datamodels.open(ref_psf_file)
        log.info(f"Using reference PSF: {ref_psf_file}")

        # Get reference filter approximate wavelength (nm) based on filter name
        ref_wavelength = get_filter_wavelength(ref_filter)

        return {
            "ref_filter": ref_filter,
            "ref_wavelength": ref_wavelength,
            "ref_model": ref_model,
            "ref_psf_model": ref_psf_model,
        }

    def _prepare_processing_order(self, library, ref_filter):
        """
        Prepare the processing order for filter models.

        Ensures the reference filter is processed first so its
        catalog is available when processing redder filters.

        Parameters
        ----------
        library : `romancal.datamodels.ModelLibrary`
            The library of models to process.

        ref_filter : str
            The reference filter name.

        Returns
        -------
        result : list of tuple
            List of (model_index, filter_name) tuples sorted so
            reference filter is first.
        """
        # Create list of model indices in processing order with
        # reference filter first, then all others.
        model_indices = []
        with library:
            for i, model in enumerate(library):
                filter_name = model.meta.instrument.optical_element
                model_indices.append((i, filter_name))
                library.shelve(model, modify=False)

        # Sort so reference filter is first
        def sort_key(index_filter_tuple):
            _, filt = index_filter_tuple
            return 0 if filt == ref_filter else 1

        model_indices.sort(key=sort_key)

        return model_indices

    def _join_filter_catalogs(self, detection_catalog, filter_catalogs):
        """
        Join filter catalogs to detection catalog in wavelength order.

        Parameters
        ----------
        detection_catalog : Table
            The detection catalog to which filter catalogs will be
            joined.

        filter_catalogs : dict
            Dictionary with (filter_name, wavelength) tuples as keys
            and catalog tables as values.

        Returns
        -------
        result : Table
            The detection catalog with all filter catalogs joined in
            wavelength order.
        """
        # Sort by wavelength (second element of tuple key)
        sorted_filter_keys = sorted(filter_catalogs.keys(), key=lambda x: x[1])

        for filter_key in sorted_filter_keys:
            cat = filter_catalogs[filter_key]
            # The outer join prevents an empty table if any
            # columns have the same name but different values
            # (e.g., repeated filter names)
            detection_catalog = join(
                detection_catalog, cat, keys="label", join_type="outer"
            )

        return detection_catalog

    def _finalize_ee_fractions(self, detection_catalog, filter_ee_fractions):
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

    def process(self, library):
        # All input MosaicImages in the ModelLibrary are assumed to
        # have the same shape and be pixel aligned.
        if isinstance(library, str):
            library = ModelLibrary(library)
        if not isinstance(library, ModelLibrary):
            msg = "library input must be a ModelLibrary object"
            raise TypeError(msg)

        with library:
            example_model = library.borrow(0)
            library.shelve(example_model, modify=False)

        # Initialize catalog model with metadata
        catalog_model = self._initialize_catalog_model(library, example_model)

        log.info("Creating ee_fractions model for first image")
        apcorr_ref = self.get_reference_file(example_model, "apcorr")
        ee_spline = get_ee_spline(example_model, apcorr_ref)

        log.info("Calculating and subtracting background")
        library = subtract_background_library(library, self.bkg_boxsize)

        # Process detection image (background, segmentation, catalog)
        detection_result = self._process_detection_image(
            library, example_model, ee_spline, catalog_model
        )

        # Check if detection failed (returns save_empty_results)
        if not isinstance(detection_result, dict):
            log.warning("Detection image processing failed.")
            return detection_result

        # Extract detection results
        segment_img = detection_result["segment_img"]
        detection_catobj = detection_result["detection_catobj"]
        detection_catalog = detection_result["detection_catalog"]
        star_kernel_fwhm = detection_result["star_kernel_fwhm"]

        # Setup reference filter for PSF matching
        ref_setup = self._setup_reference_filter(library)
        ref_filter = ref_setup["ref_filter"]
        ref_wavelength = ref_setup["ref_wavelength"]
        ref_model = ref_setup["ref_model"]
        ref_psf_model = ref_setup["ref_psf_model"]

        # Record the PSF match reference filter in metadata
        detection_catalog.meta["psf_match_reference_filter"] = ref_filter.lower()

        # Store reference filter's catalog for computing correction
        # factors
        ref_filter_catalog = None

        # Dictionary to store filter catalogs for joining in
        # wavelength order
        filter_catalogs = {}

        # List to store ee_fractions from each filter
        filter_ee_fractions = []

        # Prepare processing order (reference filter first)
        model_indices = self._prepare_processing_order(library, ref_filter)

        # Create catalogs for each input image
        time_means = []
        exposure_times = []
        with library:
            for model_index, _ in model_indices:
                model = library.borrow(model_index)
                filter_name = model.meta.instrument.optical_element

                # Create catalog for this filter
                result = create_filter_catalog(
                    model=model,
                    filter_name=filter_name,
                    ref_filter=ref_filter,
                    ref_wavelength=ref_wavelength,
                    segment_img=segment_img,
                    star_kernel_fwhm=star_kernel_fwhm,
                    detection_catobj=detection_catobj,
                    ref_model=ref_model,
                    ref_filter_catalog=ref_filter_catalog,
                    ref_psf_model=ref_psf_model,
                    fit_psf=self.fit_psf,
                    get_reference_file_func=self.get_reference_file,
                )

                # Store the reference filter catalog. This will be
                # None until we process the reference filter.
                ref_filter_catalog = result["ref_filter_catalog"]

                # Store ee_fractions for consolidation
                filter_ee_fractions.append(result["ee_fractions"])

                # Store filter catalog for later joining in wavelength
                # order
                filter_wavelength = get_filter_wavelength(filter_name)
                cat = result["catalog"]
                filter_catalogs[(filter_name, filter_wavelength)] = cat

                # Accumulate and blend image metadata
                blend_image_metadata(
                    model,
                    catalog_model,
                    time_means,
                    exposure_times,
                )

                library.shelve(model, modify=False)

        # Join all filter catalogs to detection catalog
        detection_catalog = self._join_filter_catalogs(
            detection_catalog, filter_catalogs
        )

        # Finish blending
        catalog_model.meta.coadd_info.time_mean = Time(time_means).mean()
        catalog_model.meta.coadd_info.exposure_time = np.mean(exposure_times)

        # Consolidate and sort ee_fractions
        self._finalize_ee_fractions(detection_catalog, filter_ee_fractions)

        # Put the resulting multiband catalog in the model
        catalog_model.source_catalog = detection_catalog

        return save_all_results(self, segment_img, catalog_model)
