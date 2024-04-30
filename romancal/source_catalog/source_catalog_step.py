"""
Module for the source catalog step.
"""

import numpy as np
from astropy.table import Table
from roman_datamodels import datamodels, maker_utils
from roman_datamodels.datamodels import ImageModel, MosaicModel

from romancal.source_catalog.background import RomanBackground
from romancal.source_catalog.detection import convolve_data, make_segmentation_image
from romancal.source_catalog.reference_data import ReferenceData
from romancal.source_catalog.source_catalog import RomanSourceCatalog
from romancal.stpipe import RomanStep

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
    reference_file_types = []

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
    """

    def process(self, input_model):
        with datamodels.open(input_model) as model:
            if not isinstance(model, (ImageModel, MosaicModel)):
                raise ValueError(
                    "The input model must be an ImageModel or MosaicModel."
                )

            if isinstance(model, ImageModel):
                cat_model = datamodels.SourceCatalogModel
            else:
                cat_model = datamodels.MosaicSourceCatalogModel
            source_catalog_model = maker_utils.mk_datamodel(cat_model)

            for key in source_catalog_model.meta.keys():
                try:
                    if key == "optical_element":
                        value = model.meta.instrument[key]
                    else:
                        value = model.meta[key]
                    source_catalog_model.meta[key] = value
                except KeyError:
                    pass

            aperture_ee = (self.aperture_ee1, self.aperture_ee2, self.aperture_ee3)
            refdata = ReferenceData(model, aperture_ee)
            aperture_params = refdata.aperture_params

            mask = np.isnan(model.data)
            coverage_mask = np.isnan(model.err)
            bkg = RomanBackground(
                model.data,
                box_size=self.bkg_boxsize,
                mask=mask,
                coverage_mask=coverage_mask,
            )
            model.data -= bkg.background

            convolved_data = convolve_data(
                model.data, kernel_fwhm=self.kernel_fwhm, mask=coverage_mask
            )

            segment_img = make_segmentation_image(
                convolved_data,
                snr_threshold=self.snr_threshold,
                npixels=self.npixels,
                bkg_rms=bkg.background_rms,
                deblend=self.deblend,
                mask=coverage_mask,
            )

            if segment_img is None:  # no sources found
                source_catalog_model.source_catalog = Table()
                return source_catalog_model

            ci_star_thresholds = (self.ci1_star_threshold, self.ci2_star_threshold)
            catobj = RomanSourceCatalog(
                model,
                segment_img,
                convolved_data,
                aperture_params,
                ci_star_thresholds,
                self.kernel_fwhm,
            )

            # put the resulting catalog in the model
            source_catalog_model.source_catalog = catobj.catalog

            # add back background to data so input model is unchanged
            # (in case of interactive use)
            model.data += bkg.background

            if self.save_results:
                # NOTE: the source_catalog_model is automatically saved
                #       if save_results = True

                # save the segmentation map
                segmentation_model = maker_utils.mk_datamodel(
                    datamodels.MosaicSegmentationMapModel
                )
                segmentation_model.data = segment_img.data.astype(np.uint32)
                self.save_model(segmentation_model, suffix="segm")

        return source_catalog_model
