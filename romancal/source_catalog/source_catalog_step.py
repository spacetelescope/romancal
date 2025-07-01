"""
Module for the source catalog step.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from astropy.table import join
from photutils.segmentation import SegmentationImage
from roman_datamodels import datamodels
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel
from roman_datamodels.stnode import SourceCatalog

from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import convolve_data, make_segmentation_image
from romancal.source_catalog.save_utils import save_segment_image
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(RomanStep):
    """
    Create a catalog of sources including photometry and basic shape
    measurements.

    Parameters
    -----------
    input : str, `ImageModel`, or `MosaicModel`
        Path to an ASDF file, or an `ImageModel` or `MosaicModel`.
    """

    class_alias = "source_catalog"
    reference_file_types: ClassVar = []

    spec = """
        bkg_boxsize = integer(default=1000)   # background mesh box size in pixels
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # per-pixel SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        suffix = string(default='cat')        # Default suffix for output files
        fit_psf = boolean(default=True)       # fit source PSFs for accurate astrometry?
        forced_segmentation = string(default='')  # force the use of this segmentation map
    """

    def process(self, step_input):
        if isinstance(step_input, datamodels.DataModel):
            input_model = step_input
        else:
            input_model = datamodels.open(step_input)

        if not isinstance(input_model, ImageModel | MosaicModel):
            raise ValueError("The input model must be an ImageModel or MosaicModel.")

        # get the name of the psf reference file
        self.ref_file = self.get_reference_file(input_model, "epsf")
        self.log.info("Using ePSF reference file: %s", self.ref_file)
        psf_ref_model = datamodels.open(self.ref_file)

        # Define a boolean mask for pixels to be excluded
        mask = (
            ~np.isfinite(input_model.data)
            | ~np.isfinite(input_model.err)
            | (input_model.err <= 0)
        )

        # Copy the data and error arrays to avoid modifying the input
        # model. The metadata and dq and weight arrays are not copied
        # because they are not modified in this step.
        if isinstance(input_model, ImageModel):
            model = ImageModel()
            model.meta = input_model.meta
            model.data = input_model.data.copy()
            model.err = input_model.err.copy()
            model.dq = input_model.dq

            # Create a DQ mask for pixels to be excluded; currently all
            # pixels with any DQ flag are excluded from the source catalog
            # except for those in ignored_dq_flags.
            # TODO: revisit these flags when CRDS reference files are updated
            ignored_dq_flags = pixel.NO_LIN_CORR
            dq_mask = np.any(model.dq[..., None] & ~ignored_dq_flags, axis=-1)

            # TODO: to set the mask to True for *only* dq_flags use:
            # dq_mask = np.any(model.dq[..., None] & dq_flags, axis=-1)
            mask |= dq_mask
        elif isinstance(input_model, MosaicModel):
            model = MosaicModel()
            model.meta = input_model.meta
            model.data = input_model.data.copy()
            model.err = input_model.err.copy()
            model.weight = input_model.weight

        # Initialize the source catalog model, copying the metadata
        # from the input model
        if isinstance(model, ImageModel):
            if self.forced_segmentation:
                cat_model_cls = datamodels.ForcedImageSourceCatalogModel
            else:
                cat_model_cls = datamodels.ImageSourceCatalogModel
        else:
            if self.forced_segmentation:
                cat_model_cls = datamodels.ForcedMosaicSourceCatalogModel
            else:
                cat_model_cls = datamodels.MosaicSourceCatalogModel
        cat_model = cat_model_cls.create_minimal({"meta": model.meta})
        if "instrument" in model.meta:
            cat_model.meta.optical_element = model.meta.instrument.optical_element
        if self.forced_segmentation:
            cat_model.meta["forced_segmentation"] = self.forced_segmentation

        # Return an empty segmentation image and catalog table if all
        # pixels are masked
        if np.all(mask):
            self.log.warning("Cannot create source catalog. All pixels are masked.")
            segment_img = np.zeros(model.data.shape, dtype=np.uint32)
            cat_model.source_catalog = cat_model.create_empty_catalog()
            return self.save_all_results(input_model, segment_img, cat_model)

        self.log.info("Calculating and subtracting background")
        bkg = RomanBackground(
            model.data,
            box_size=self.bkg_boxsize,
            coverage_mask=mask,
        )
        model.data -= bkg.background

        self.log.info("Creating detection image")
        detection_image = convolve_data(
            model.data, kernel_fwhm=self.kernel_fwhm, mask=mask
        )

        self.log.info("Detecting sources")
        if not self.forced_segmentation:
            segment_img = make_segmentation_image(
                detection_image,
                snr_threshold=self.snr_threshold,
                npixels=self.npixels,
                bkg_rms=bkg.background_rms,
                deblend=self.deblend,
                mask=mask,
            )
            if segment_img is not None:
                segment_img.detection_image = detection_image
        else:
            forced_segmodel = datamodels.open(self.forced_segmentation)
            # forced_segmodel.data is asdf.tags.core.ndarray.NDArrayType
            forced_segimg = forced_segmodel.data[...]

            # Remove fully masked segments
            unmasked_sources = np.unique(forced_segimg * (mask == 0))
            fully_masked_sources = set(np.unique(forced_segimg)) - set(unmasked_sources)
            forced_segimg_mask = np.isin(
                forced_segimg, np.array(list(fully_masked_sources))
            )
            forced_segimg[forced_segimg_mask] = 0
            segment_img = SegmentationImage(forced_segimg)

        # Return an empty segmentation image and catalog table if no
        # sources are detected
        if segment_img is None:
            self.log.warning("Cannot create source catalog. No sources were detected.")
            segment_img = np.zeros(model.data.shape, dtype=np.uint32)
            cat_model.source_catalog = cat_model.create_empty_catalog()
            return self.save_all_results(input_model, segment_img, cat_model)

        self.log.info("Creating source catalog")
        cat_type = "prompt" if not self.forced_segmentation else "forced_det"
        fit_psf = self.fit_psf & (not self.forced_segmentation)  # skip when forced
        catobj = RomanSourceCatalog(
            model,
            segment_img,
            detection_image,
            self.kernel_fwhm,
            fit_psf=fit_psf,
            mask=mask,
            psf_ref_model=psf_ref_model,
            cat_type=cat_type,
        )
        cat = catobj.catalog

        if self.forced_segmentation:
            # TODO: improve this so that the moment-based properties are
            # not recomputed from the forced_detection_image
            forced_detection_image = forced_segmodel.detection_image
            segment_img.detection_image = forced_detection_image
            forced_catobj = RomanSourceCatalog(
                model,
                segment_img,
                forced_detection_image,
                self.kernel_fwhm,
                fit_psf=self.fit_psf,
                mask=mask,
                psf_ref_model=psf_ref_model,
                cat_type="forced_full",
            )

            # We have two catalogs, both using the same segmentation
            # image. We want:
            # - the original shape parameters computed from
            #   the forced detection image.  These are needed to
            #   describe where we have computed the forced photometry.
            #   These keep their original names to match up with the deep
            #   catalog used for forcing.
            # - the newly measured fluxes and flags and sharpness
            #   / roundness from the direct image; these give the new fluxes
            #   at these locations
            #   These gain a forced_ prefix.
            # - the shapes measured from the new detection image.  These
            #   seem to me to have less value but are explicitly called out in
            #   a requirement, and it's not crazy to compute new centroids and
            #   moments.
            #   These gain a forced_prefix.
            # At the end of the day you get a whole new catalog with the forced_
            # prefix, plus some shape parameters that duplicate values in the
            # original catalog used for forcing.

            # merge the two forced catalogs
            forced_cat = forced_catobj.catalog
            forced_cat.meta = None  # redundant with cat.meta
            cat = join(forced_cat, cat, keys="label", join_type="outer")

        # Put the resulting catalog table in the catalog model
        cat_model.source_catalog = cat

        return self.save_all_results(input_model, segment_img, cat_model)

    def save_all_results(self, model, segment_img, cat_model):
        """
        Return and save the results of the source catalog step.

        The segmentation image is always saved.

        If ``save_results`` is True, then either the input model with
        updated metadata or the source catalog model is saved/returned
        depending on the value of ``return_updated_model``. If
        ``return_updated_model`` is True, then the input model is
        saved/returned with updated metadata. If False, then the source
        catalog model is saved/returned.

        The source catalog is saved as a parquet file.

        Parameters
        ----------
        model : `ImageModel` or `MosaicModel`
            The input model to the source catalog step.

        segment_img : `photutils.segmentation.SegmentationImage`
            The segmentation image created by the source catalog step.

        cat_model : `datamodels.ImageSourceCatalogModel`, `datamodels.MosaicSourceCatalogModel`, `datamodels.ForcedImageSourceCatalogModel`, or `datamodels.ForcedMosaicSourceCatalogModel`
            The source catalog model created by the source catalog step.

        Returns
        -------
        result : `ImageModel`, `MosaicModel`, or `SourceCatalog`
            The source catalog model or the input model with updated
            metadata.
        """
        if isinstance(segment_img, np.ndarray):
            # convert the segmentation image to a SegmentationImage
            segment_img = SegmentationImage(segment_img)

        # define the output filename
        output_filename = (
            self.output_file
            if self.output_file is not None
            else cat_model.meta.filename
        )

        # always save the segmentation image
        save_segment_image(self, segment_img, cat_model, output_filename)

        # Update the source catalog filename metadata
        self.output_ext = "parquet"
        output_catalog_name = self.make_output_path(
            basepath=output_filename, suffix="cat"
        )
        self.output_ext = "asdf"

        # Always save the source catalog, but don't save it twice.
        # If save_results=False or return_update_model=True, we need to
        # explicitly save it.
        return_updated_model = getattr(self, "return_updated_model", False)
        if not self.save_results or return_updated_model:
            self.output_ext = "parquet"
            self.save_model(
                cat_model,
                output_file=output_filename,
                suffix="cat",
                force=True,
            )
            self.output_ext = "asdf"

        # Return the source catalog object or the input model. If the
        # input model is an ImageModel, the metadata is updated with the
        # source catalog filename.
        if getattr(self, "return_updated_model", False):
            # define the catalog filename; self.save_model will
            # determine whether to use a fully qualified path

            # set the suffix to something else to prevent the step from
            # overwriting the source catalog file with a datamodel
            self.suffix = "sourcecatalog"

            # update the input datamodel with the tweakreg catalog name
            if isinstance(model, ImageModel):
                model.meta.source_catalog = SourceCatalog()
                model.meta.source_catalog.tweakreg_catalog_name = output_catalog_name
                model.meta.cal_step.source_catalog = "COMPLETE"

            result = model
        else:
            self.output_ext = "parquet"
            result = cat_model

        # validate the result to flush out any lazy-loaded contents
        result.validate()
        return result
