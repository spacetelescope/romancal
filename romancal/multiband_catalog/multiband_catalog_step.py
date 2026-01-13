"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING

from roman_datamodels import datamodels

from romancal.datamodels.fileio import open_dataset
from romancal.multiband_catalog.multiband_catalog import (
    make_source_injected_library,
    match_recovered_sources,
    multiband_catalog,
)
from romancal.source_catalog.save_utils import save_all_results, save_empty_results
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
        snr_threshold = float(default=5.0)    # per-pixel SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        suffix = string(default='cat')        # Default suffix for output files
        fit_psf = boolean(default=True)       # fit source PSFs for accurate astrometry?
        inject_sources = boolean(default=False) # Inject sources into images
        save_debug_info = boolean(default=False)
                                   # Include image data and other data for testing
    """

    def process(self, dataset):
        # All input MosaicImages in the ModelLibrary are assumed to have
        # the same shape and be pixel aligned.
        library = open_dataset(
            dataset, update_version=self.update_version, as_library=True
        )

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

        # Set up source injection library and injection catalog
        if self.inject_sources:
            si_library, si_cat = make_source_injected_library(library)

        # Create catalog of library images
        segment_img, cat_model, msg = multiband_catalog(
            self, library, example_model, cat_model, ee_spline
        )

        # The results are empty
        if msg is not None:
            return save_empty_results(self, segment_img, cat_model, msg=msg)

        # Source Injection
        if self.inject_sources:
            with si_library:
                si_example_model = si_library.borrow(0)
                si_library.shelve(si_example_model, modify=False)

            si_ee_spline = get_ee_spline(si_example_model, apcorr_ref)

            # Create catalog of source injected images
            si_segment_img, si_cat_model, _ = multiband_catalog(
                self,
                si_library,
                si_example_model,
                copy.deepcopy(cat_model),
                si_ee_spline,
            )

            # Match sources
            recovered_sources = match_recovered_sources(
                cat_model.source_catalog, si_cat, si_cat_model.source_catalog
            )

            # Put the source injected multiband catalog in the model
            cat_model.source_injection_catalog = si_cat_model.source_catalog
            segment_img.injected_sources = si_cat
            segment_img.recovered_sources = recovered_sources

            if self.save_debug_info:
                segment_img.si_segment_img = si_segment_img
                segment_img.si_detection_image = si_segment_img.detection_image

        return save_all_results(
            self, segment_img, cat_model, save_debug_info=self.save_debug_info
        )
