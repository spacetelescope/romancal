"""
Module for the multiband source catalog step.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from astropy.table import Table, join
from roman_datamodels import datamodels

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog.background import subtract_background_library
from romancal.multiband_catalog.detection_image import make_detection_image
from romancal.multiband_catalog.utils import add_filter_to_colnames
from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import make_segmentation_image
from romancal.source_catalog.save_utils import save_segment_image
from romancal.source_catalog.source_catalog import RomanSourceCatalog
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
    reference_file_types: ClassVar = []

    spec = """
        bkg_boxsize = integer(default=100)   # background mesh box size in pixels
        kernel_fwhms = float_list(default=None)  # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # per-pixel SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        suffix = string(default='cat')        # Default suffix for output files
        fit_psf = boolean(default=True)       # fit source PSFs for accurate astrometry?
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

        source_catalog_model = datamodels.MultibandSourceCatalogModel.create_minimal(
            {"meta": example_model.meta}
        )
        if "instrument" in example_model.meta:
            source_catalog_model.meta.optical_element = (
                example_model.meta.instrument.optical_element
            )

        try:
            source_catalog_model.meta.filename = library.asn["products"][0]["name"]
        except (AttributeError, KeyError):
            source_catalog_model.meta.filename = "multiband_catalog"
        if self.output_file is None:
            self.output_file = source_catalog_model.meta.filename

        # TODO: sensible defaults
        # TODO: redefine in terms of intrinsic FWHM
        if self.kernel_fwhms is None:
            self.kernel_fwhms = [2.0, 20.0]

        library = subtract_background_library(library, self.bkg_boxsize)

        # TODO: where to save the det_img and det_err?
        det_img, det_err = make_detection_image(library, self.kernel_fwhms)

        # estimate background rms from detection image to calculate a
        # threshold for source detection
        mask = ~np.isfinite(det_img) | ~np.isfinite(det_err) | (det_err <= 0)
        bkg = RomanBackground(
            det_img,
            box_size=self.bkg_boxsize,
            coverage_mask=mask,
        )
        bkg_rms = bkg.background_rms

        segment_img = make_segmentation_image(
            det_img,
            snr_threshold=self.snr_threshold,
            npixels=self.npixels,
            bkg_rms=bkg_rms,
            deblend=self.deblend,
            mask=mask,
        )
        segment_img.detection_image = det_img

        if segment_img is None:  # no sources found
            source_catalog_model.source_catalog = Table()
        else:
            # this is needed for the DAOStarFinder sharpness and roundness
            # properties
            # TODO: measure on a secondary detection image with minimal
            # smoothing; same for basic shape parameters
            star_kernel_fwhm = np.min(self.kernel_fwhms)  # ??

            det_model = datamodels.MosaicModel()
            det_model.data = det_img
            det_model.err = det_err

            # TODO: this is a temporary solution to get model attributes
            # currently needed in RomanSourceCatalog
            det_model.weight = example_model.weight
            det_model.meta = example_model.meta

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
            )
            # need to generate the catalog before we pass the det_catobj
            # to the RomanSourceCatalog constructor
            det_cat = det_catobj.catalog

            # loop over each image
            with library:
                for model in library:
                    mask = (
                        ~np.isfinite(model.data)
                        | ~np.isfinite(model.err)
                        | (model.err <= 0)
                    )

                    if self.fit_psf:
                        filter_name = model.meta.basic.optical_element  # L3
                        log.info(f"Creating catalog for {filter_name} image")
                        ref_file = self.get_reference_file(model, "epsf")
                        self.log.info("Using ePSF reference file: %s", ref_file)
                        psf_ref_model = datamodels.open(ref_file)
                    else:
                        psf_ref_model = None

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
                    )

                    filter_name = model.meta.basic.optical_element
                    cat = add_filter_to_colnames(catobj.catalog, filter_name)
                    # TODO: what metadata do we want to keep, if any,
                    # from the filter catalogs?
                    cat.meta = None  # temporary
                    # outer join prevents empty table if any columns have
                    # the same name but different values (e.g., repeated
                    # filter names)
                    det_cat = join(det_cat, cat, keys="label", join_type="outer")
                    library.shelve(model, modify=False)

            # put the resulting catalog in the model
            source_catalog_model.source_catalog = det_cat

        # always save the segmentation image
        output_filename = (
            self.output_file
            if self.output_file is not None
            else source_catalog_model.meta.filename
        )
        save_segment_image(self, segment_img, source_catalog_model, output_filename)

        self.output_ext = "parquet"

        return source_catalog_model
