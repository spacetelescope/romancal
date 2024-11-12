"""
Module for the multiband source catalog step.
"""

import logging

import numpy as np
from astropy.table import Table, join
from roman_datamodels import datamodels, maker_utils

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog.background import subtract_background_library
from romancal.multiband_catalog.detection_image import make_detection_image
from romancal.multiband_catalog.utils import update_colnames
from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import make_segmentation_image
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.stpipe import RomanStep

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
    reference_file_types = []

    spec = """
        bkg_boxsize = integer(default=1000)   # background mesh box size in pixels
        kernel_fwhms = float_list(default=None)  # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # per-pixel SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        aperture_ee1 = integer(default=30)    # aperture encircled energy 1
        aperture_ee2 = integer(default=50)    # aperture encircled energy 2
        aperture_ee3 = integer(default=70)    # aperture encircled energy 3
        ci1_star_threshold = float(default=2.0)  # CI 1 star threshold
        ci2_star_threshold = float(default=1.8)  # CI 2 star threshold
        suffix = string(default='cat')        # Default suffix for output files
        fit_psf = boolean(default=True)      # fit source PSFs for accurate astrometry?
    """

    def process(self, library):
        # All input MosaicImages in the ModelLibrary are assumed to have
        # the same shape and be pixel aligned.
        if isinstance(library, str):
            library = ModelLibrary(library)
        if not isinstance(library, ModelLibrary):
            raise TypeError("library input must be a ModelLibrary object")

        # TODO: sensible defaults
        if self.kernel_fwhms is None:
            self.kernel_fwhms = [2.0, 20.0]

        # NOTE: I'm assuming here that all of the input images have
        # have the same shape and are pixel aligned.

        cat_model = datamodels.MosaicSourceCatalogModel
        source_catalog_model = maker_utils.mk_datamodel(cat_model)
        source_catalog_model.meta.filename = library.asn["products"][0]["name"]

        # define the aperture parameters for the source catalog
        # based on the input encircled energy fractions
        # TODO: where will these values be defined?
        # source_catalog will ultimately get these filter-dependent
        # values from a reference file based on EE values;
        # do we want filter-dependent aperture parameters for the
        # multiband catalog?
        aperture_ee = (self.aperture_ee1, self.aperture_ee2, self.aperture_ee3)
        aperture_params = {
            "aperture_radii": np.array((1.0, 2.0, 3.0)),
            "aperture_corrections": np.array((1.0, 1.0, 1.0)),
            "aperture_ee": aperture_ee,
            "bkg_aperture_inner_radius": 5.0,
            "bkg_aperture_outer_radius": 10.0,
        }

        library = subtract_background_library(library, self.bkg_boxsize)

        # TODO: do we want to save the det_img and det_err?
        det_img, det_err = make_detection_image(library, self.kernel_fwhms)

        # estimate background rms from detection image to calculate a
        # threshold for source detection
        mask = np.isnan(det_img)
        coverage_mask = np.isnan(det_err)
        bkg = RomanBackground(
            det_img,
            box_size=self.bkg_boxsize,
            mask=mask,
            coverage_mask=coverage_mask,
        )
        bkg_rms = bkg.background_rms

        segment_img = make_segmentation_image(
            det_img - np.nanmedian(det_img),
            snr_threshold=self.snr_threshold,
            npixels=self.npixels,
            bkg_rms=bkg_rms,
            deblend=self.deblend,
            mask=coverage_mask,
        )

        if segment_img is None:  # no sources found
            source_catalog_model.source_catalog = Table()
        else:
            ci_star_thresholds = (
                self.ci1_star_threshold,
                self.ci2_star_threshold,
            )

            # this is needed for the DAOFind sharpness and roundness
            # properties; are these needed for the Roman source catalog?
            star_kernel_fwhm = np.min(self.kernel_fwhms)  # ??

            det_model = maker_utils.mk_datamodel(datamodels.MosaicModel)
            det_model.data = det_img
            det_model.err = det_err

            # TODO: this is a temporary solution to get model attributes
            # currently needed in RomanSourceCatalog
            with library:
                model = library.borrow(0)
                det_model.weight = model.weight
                det_model.meta = model.meta
                library.shelve(model, modify=False)

            det_catobj = RomanSourceCatalog(
                det_model,
                segment_img,
                det_img,
                aperture_params,
                ci_star_thresholds,
                star_kernel_fwhm,
                self.fit_psf,
                detection_cat=None,
            )
            det_cat = update_colnames(det_catobj.catalog, prefix="det_")

            # loop over each image
            with library:
                for model in library:
                    catobj = RomanSourceCatalog(
                        model,
                        segment_img,
                        None,
                        aperture_params,
                        ci_star_thresholds,
                        star_kernel_fwhm,
                        self.fit_psf,
                        detection_cat=det_catobj,
                    )

                    filter_name = model.meta.basic.optical_element
                    prefix = f"{filter_name}_"
                    cat = update_colnames(catobj.catalog, prefix=prefix)
                    cat.meta = None  # temporary
                    # outer join prevents empty table if any columns have
                    # the same name but different values (e.g., repeated
                    # filter names)
                    det_cat = join(det_cat, cat, keys="label", join_type="outer")
                    library.shelve(model, modify=False)

            # put the resulting catalog in the model
            source_catalog_model.source_catalog = det_cat

        # save the segmentation image and multiband catalog
        # TODO: I noticed that the catalog is saved twice;
        # once here and once when the step returns
        self.save_base_results(segment_img, source_catalog_model)

        return source_catalog_model

    def save_base_results(self, segment_img, source_catalog_model):
        # save the segmentation map and source catalog
        output_filename = (
            self.output_file
            if self.output_file is not None
            else source_catalog_model.meta.filename
        )

        seg_model = datamodels.MosaicSegmentationMapModel
        segmentation_model = maker_utils.mk_datamodel(seg_model)
        for key in segmentation_model.meta.keys():
            segmentation_model.meta[key] = source_catalog_model.meta[key]

        if segment_img is not None:
            segmentation_model.data = segment_img.data.astype(np.uint32)
            self.save_model(
                segmentation_model,
                output_file=output_filename,
                suffix="segm",
                force=True,
            )
        # save the source catalog
        self.save_model(
            source_catalog_model,
            output_file=output_filename,
            suffix="cat",
            force=True,
        )
