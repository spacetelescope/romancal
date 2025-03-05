"""
Module for the source catalog step.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from astropy.table import Table, join
from photutils.segmentation import SegmentationImage
from roman_datamodels import datamodels, maker_utils
from roman_datamodels.datamodels import ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel
from roman_datamodels.maker_utils import mk_datamodel

from romancal.multiband_catalog import utils
from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import convolve_data, make_segmentation_image
from romancal.source_catalog.reference_data import ReferenceData
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
        aperture_ee1 = integer(default=30)    # aperture encircled energy 1
        aperture_ee2 = integer(default=50)    # aperture encircled energy 2
        aperture_ee3 = integer(default=70)    # aperture encircled energy 3
        ci1_star_threshold = float(default=2.0)  # CI 1 star threshold
        ci2_star_threshold = float(default=1.8)  # CI 2 star threshold
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

        # Define a boolean mask for pixels to be excluded
        mask = (
            ~np.isfinite(input_model.data)
            | ~np.isfinite(input_model.err)
            | (input_model.err <= 0)
        )

        # Copy the data and error arrays to avoid modifying the input
        # model. We use mk_datamodel to copy *only* the data and err
        # arrays. The metadata and dq and weight arrays are not copied
        # because they are not modified in this step. The other model
        # arrays (e.g., var_rnoise) are not currently used by this step.
        if isinstance(input_model, ImageModel):
            model = mk_datamodel(
                ImageModel,
                meta=input_model.meta,
                shape=(0, 0),
                data=input_model.data.copy(),
                err=input_model.err.copy(),
                dq=input_model.dq,
            )

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
            model = mk_datamodel(
                MosaicModel,
                meta=input_model.meta,
                shape=(0, 0),
                data=input_model.data.copy(),
                err=input_model.err.copy(),
                weight=input_model.weight,
            )

        if isinstance(model, ImageModel):
            cat_model = datamodels.ImageSourceCatalogModel
        else:
            cat_model = datamodels.MosaicSourceCatalogModel
        source_catalog_model = maker_utils.mk_datamodel(cat_model)

        for key in source_catalog_model.meta.keys():
            value = (
                model.meta.instrument[key]
                if key == "optical_element"
                else model.meta[key]
            )
            source_catalog_model.meta[key] = value

        if self.forced_segmentation != "":
            source_catalog_model.meta["forced_segmentation"] = self.forced_segmentation
            forced = True
        else:
            forced = False

        aperture_ee = (self.aperture_ee1, self.aperture_ee2, self.aperture_ee3)
        refdata = ReferenceData(model, aperture_ee)
        aperture_params = refdata.aperture_params

        bkg = RomanBackground(
            model.data,
            box_size=self.bkg_boxsize,
            coverage_mask=mask,
        )
        model.data -= bkg.background

        detection_image = convolve_data(
            model.data, kernel_fwhm=self.kernel_fwhm, mask=mask
        )

        if not forced:
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
            segment_img = SegmentationImage(forced_segmodel.data[...])

        if segment_img is None:  # no sources found
            cat = Table()
        else:
            ci_star_thresholds = (
                self.ci1_star_threshold,
                self.ci2_star_threshold,
            )
            catobj = RomanSourceCatalog(
                model,
                segment_img,
                detection_image,
                aperture_params,
                ci_star_thresholds,
                self.kernel_fwhm,
                self.fit_psf & (not forced),
                # don't need to do PSF photometry here when forcing; happens later
                mask=mask,
            )
            cat = catobj.catalog

        if forced:
            forced_detection_image = forced_segmodel.detection_image[...]
            segment_img.detection_image = forced_detection_image
            forcedcatobj = RomanSourceCatalog(
                model,
                segment_img,
                forced_detection_image,
                aperture_params,
                ci_star_thresholds,
                self.kernel_fwhm,
                self.fit_psf,
                mask=mask,
            )
            # we have two catalogs, both using a specified set of
            # pre-specified segments.  We want:
            # - keep the original shape parameters computed from
            #   the forced detection image.  These are needed to
            #   describe where we have computed the forced photometry.
            #   These keep their original names to match up with the deep
            #   catalog used for forcing.
            # - keep the newly measured fluxes and flags and sharpness
            #   / roundness from the direct image; these give the new fluxes
            #   at these locations
            #   These gain a forced_ prefix.
            # - keep the shapes measured from the new detection image.  These
            #   seem to me to have less value but are explicitly called out in
            #   a requirement, and it's not crazy to compute new centroids and
            #   moments.
            #   These gain a forced_prefix.
            # - discard fluxes we measure with the new shapes from the new
            #   detection image.
            # At the end of the day you get a whole new catalog with the forced_
            # prefix, plus some shape parameters that duplicate values in the
            # original catalog used for forcing.
            forcedcat = forcedcatobj.catalog
            utils.prefix_colnames(
                forcedcat, "forced_", colnames=utils.get_direct_image_columns(forcedcat)
            )
            utils.prefix_colnames(
                cat, "forced_", colnames=utils.get_detection_image_columns(cat)
            )
            forcedcat.remove_columns([x for x in forcedcat.colnames if "_bbox_" in x])
            cat.remove_columns(utils.get_direct_image_columns(cat))
            forcedcat.meta = None  # redundant with cat.meta
            cat = join(forcedcat, cat, keys="label", join_type="outer")
            colnames = [x for x in cat.colnames if not x.startswith("forced_")] + [
                x for x in cat.colnames if x.startswith("forced_")
            ]
            cat = cat[colnames]

        # put the resulting catalog in the model
        source_catalog_model.source_catalog = cat

        # always save the segmentation image and source catalog
        self.save_base_results(segment_img, source_catalog_model)

        # Return the source catalog object or the input model. If the
        # input model is an ImageModel, the metadata is updated with the
        # source catalog filename.
        if getattr(self, "return_updated_model", False):
            # define the catalog filename; self.save_model will
            # determine whether to use a fully qualified path
            output_catalog_name = self.make_output_path(
                basepath=model.meta.filename, suffix="cat"
            )

            # set the suffix to something else to prevent the step from
            # overwriting the source catalog file with a datamodel
            self.suffix = "sourcecatalog"

            if isinstance(input_model, ImageModel):
                update_metadata(input_model, output_catalog_name)

            result = input_model
        else:
            result = source_catalog_model

        return result

    def save_base_results(self, segment_img, source_catalog_model):
        # save the segmentation map and source catalog
        output_filename = (
            self.output_file
            if self.output_file is not None
            else source_catalog_model.meta.filename
        )

        if isinstance(source_catalog_model, datamodels.ImageSourceCatalogModel):
            seg_model = datamodels.SegmentationMapModel
        else:
            seg_model = datamodels.MosaicSegmentationMapModel

        segmentation_model = maker_utils.mk_datamodel(seg_model)
        for key in segmentation_model.meta.keys():
            segmentation_model.meta[key] = source_catalog_model.meta[key]

        if segment_img is not None:
            segmentation_model.data = segment_img.data.astype(np.uint32)
            segmentation_model["detection_image"] = segment_img.detection_image
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


def update_metadata(model, output_catalog_name):
    # update datamodel to point to the source catalog file destination
    model.meta["source_catalog"] = maker_utils.mk_source_catalog(
        tweakreg_catalog_name=output_catalog_name
    )
    model.meta.cal_step.source_catalog = "COMPLETE"
