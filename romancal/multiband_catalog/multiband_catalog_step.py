"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from roman_datamodels import datamodels

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog.multiband_catalog import multiband_catalog
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
        snr_threshold = float(default=3.0)    # per-pixel SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        suffix = string(default='cat')        # Default suffix for output files
        fit_psf = boolean(default=True)       # fit source PSFs for accurate astrometry?
        inject_sources = boolean(default=False) # Inject sources into images
        test_data = boolean(default=False)
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

        # Set ups source injection files and library
        if self.inject_sources:
            # Obtain exposure times and filters
            # This code assumes all filters have been coadded already,
            # and thus there is one image per filter
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
                    si_exptimes[si_filter_name] = \
                        float(si_model.meta.coadd_info.exposure_time)
                    si_filters.append(si_filter_name)

                    # Poisson variance required for source injection
                    if "var_poisson" not in si_model:
                        si_model.var_poisson = si_model.err**2

                    # Set parameters for source injection
                    # This only needs to be done once per library
                    if si_cen is None:
                        # Create source grid points
                        si_x_pos, si_y_pos = injection.make_source_grid(si_model,
                            yxmax=si_model.data.shape, yxoffset=(50, 50),
                            yxgrid=(20, 20))

                        si_cen = SkyCoord(ra=si_model.meta.wcsinfo.ra_ref * u.deg,
                                            dec=si_model.meta.wcsinfo.dec_ref * u.deg,)

                        # Convert to RA & Dec
                        wcsobj = si_model.meta.wcs
                        si_ra, si_dec = wcsobj.pixel_to_world_values(np.array(si_x_pos),
                            np.array(si_y_pos))

                        # Generate cosmos-like catalog
                        si_cat = injection.make_cosmoslike_catalog(
                            cen=si_cen, ra=si_ra, dec=si_dec, exptimes=si_exptimes,
                        )

                    # Inject sources into the detection image
                    si_model = injection.inject_sources(si_model, si_cat)

                    # Add model to list for new library
                    si_model_lst.append(si_model)

            # Create library of source injected models
            si_library = ModelLibrary(si_model_lst)

        # Create catalog of library images
        segment_img, cat_model, msg = multiband_catalog(self,
            library, example_model, cat_model, ee_spline)

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
            si_segment_img, si_cat_model, _ = multiband_catalog(self,
                si_library, si_example_model, copy.deepcopy(cat_model), si_ee_spline)

            # Put the source injected multiband catalog in the model
            cat_model.source_injection_catalog = si_cat_model.source_catalog
            segment_img.injected_sources = si_cat.as_array()

            if self.test_data:
                segment_img.si_segment_img = si_segment_img
                segment_img.si_detection_image = si_segment_img.detection_image

        return save_all_results(self, segment_img, cat_model, test_data=self.test_data)
