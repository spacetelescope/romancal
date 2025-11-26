"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog.multiband_catalog import (
    initialize_catalog_model,
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

        # Initialize the source catalog model
        cat_model = initialize_catalog_model(library, example_model)

        log.info("Creating ee_fractions model for first image")
        apcorr_ref = self.get_reference_file(example_model, "apcorr")
        ee_spline = get_ee_spline(example_model, apcorr_ref)

        # Set up source injection library and injection catalog
        if self.inject_sources:
            si_library, si_cat = make_source_injected_library(library)

        # Create the multiband catalog
        *results, msg = multiband_catalog(
            self, library, example_model, cat_model, ee_spline
        )

        # Save empty results if there was an error
        if msg is None:
            segment_img, cat_model = results
        else:
            return save_empty_results(self, *results, msg=msg)

        # Source Injection
        if self.inject_sources:
            with si_library:
                si_example_model = si_library.borrow(0)
                si_library.shelve(si_example_model, modify=False)

            # Create catalog of source injected images
            si_segment_img, si_cat_model, _ = multiband_catalog(
                self,
                si_library,
                si_example_model,
                copy.deepcopy(cat_model),
                ee_spline,
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
