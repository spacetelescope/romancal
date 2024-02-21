"""
Module for the source catalog step.
"""

import numpy as np
from astropy.table import Table
from crds.core.exceptions import CrdsLookupError
from roman_datamodels import datamodels, maker_utils

from romancal.stpipe import RomanStep

from .detection import RomanBackground, RomanSourceFinder, convolve_data
from .reference_data import ReferenceData
from .source_catalog import RomanSourceCatalog

__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(RomanStep):
    """
    Create a final catalog of source photometry and morphologies.

    Parameters
    -----------
    input : str or `ImageModel`
        Path to an ASDF file, or an `ImageModel`.
    """

    class_alias = "source_catalog"

    spec = """
        bkg_boxsize = integer(default=1000)   # background mesh box size in pixels
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        aperture_ee1 = integer(default=30)    # aperture encircled energy 1
        aperture_ee2 = integer(default=50)    # aperture encircled energy 2
        aperture_ee3 = integer(default=70)    # aperture encircled energy 3
        ci1_star_threshold = float(default=2.0)  # CI 1 star threshold
        ci2_star_threshold = float(default=1.8)  # CI 2 star threshold
        suffix = string(default='cat')        # Default suffix for output files
    """

    def _get_reffile_paths(self, model):
        filepaths = []
        for reffile_type in self.reference_file_types:
            try:
                filepath = self.get_reference_file(model, reffile_type)
                self.log.info(
                    f"Using {reffile_type.upper()} reference file: " f"{filepath}"
                )
            except CrdsLookupError as err:
                msg = f"{err} Source catalog will not be created."
                self.log.warning(msg)
                return None

            filepaths.append(filepath)
        return filepaths

    def process(self, input_model):
        with datamodels.open(input_model) as model:
            reffile_paths = self._get_reffile_paths(model)
            aperture_ee = (self.aperture_ee1, self.aperture_ee2, self.aperture_ee3)

            try:
                refdata = ReferenceData(model, reffile_paths, aperture_ee)
                aperture_params = refdata.aperture_params

                # TODO (@bmorris3): replace or remove
                abvega_offset = refdata.abvega_offset
            except RuntimeError as err:
                msg = f"{err} Source catalog will not be created."
                self.log.warning(msg)
                return None

            coverage_mask = np.isnan(model.err)  # | (model.wht == 0)
            bkg = RomanBackground(
                model.data, box_size=self.bkg_boxsize, coverage_mask=coverage_mask
            )
            model.data -= bkg.background

            threshold = self.snr_threshold * bkg.background_rms
            finder = RomanSourceFinder(threshold, self.npixels, deblend=self.deblend)

            convolved_data = convolve_data(
                model.data, self.kernel_fwhm, mask=coverage_mask
            )
            segment_img = finder(convolved_data, mask=coverage_mask)

            source_catalog_model = maker_utils.mk_datamodel(
                datamodels.SourceCatalogModel
            )

            if segment_img is None:
                source_catalog_model.source_catalog = Table()
                return source_catalog_model

            ci_star_thresholds = (self.ci1_star_threshold, self.ci2_star_threshold)
            catobj = RomanSourceCatalog(
                model,
                segment_img,
                convolved_data,
                self.kernel_fwhm,
                aperture_params,
                abvega_offset,
                ci_star_thresholds,
            )
            catalog = catobj.catalog

            # add back background to data so input model is unchanged
            model.data += bkg.background

            if self.save_results:
                cat_filepath = self.make_output_path(ext=".ecsv")
                catalog.write(cat_filepath, format="ascii.ecsv", overwrite=True)
                self.log.info(f"Wrote source catalog: {cat_filepath}")

                segmentation_model = maker_utils.mk_datamodel(
                    datamodels.SegmentationMapModel
                )
                self.save_model(segmentation_model, suffix="segm")

                # TODO (@bmorris3): replace or remove
                # model.meta.segmentation_map = segmentation_model.meta.filename

                self.log.info(
                    "Wrote segmentation map: " f"{segmentation_model.meta.filename}"
                )

        # put the resulting catalog in the model:
        source_catalog_model.source_catalog = catalog

        return source_catalog_model
