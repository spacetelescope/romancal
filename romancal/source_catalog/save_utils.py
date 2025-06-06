import logging

import numpy as np
from roman_datamodels.datamodels import (
    ForcedImageSourceCatalogModel,
    ImageSourceCatalogModel,
    MosaicSegmentationMapModel,
    SegmentationMapModel,
)


def save_segment_image(self, segment_img, source_catalog_model, output_filename):
    """
    Save the segmentation image to an ASDF file.

    Parameters
    ----------
    self : `romancal.source_catalog.SourceCatalogStep` or `romancal.multiband_catalog.MultiBandCatalogStep`
        The step instance that is calling this function.

    segment_img : `photutils.segmentation.SegmentationImage`
        The segmentation image to save.

    source_catalog_model : `roman_datamodels.datamodels.ImageSourceCatalogModel` or `roman_datamodels.datamodels.MosaicSourceCatalogModel` or `roman_datamodels.datamodels.ForcedImageSourceCatalogModel` or `roman_datamodels.datamodels.ForcedMosaicSourceCatalogModel`
        The source catalog model to use for metadata.

    output_filename : str
        The output file name.
    """
    if segment_img is None:
        logging.warning("No segmentation image to save.")
        return

    # Define the segmentation model based on the source catalog model
    # type
    if isinstance(
        source_catalog_model, (ImageSourceCatalogModel | ForcedImageSourceCatalogModel)
    ):
        segm_model = SegmentationMapModel
    else:
        segm_model = MosaicSegmentationMapModel
    segmentation_model = segm_model.create_minimal({"meta": source_catalog_model.meta})

    # Set the data and detection image
    segmentation_model.data = segment_img.data.astype(np.uint32)
    segmentation_model["detection_image"] = segment_img.detection_image

    # Save the segmentation image to the output file
    self.output_ext = "asdf"
    self.save_model(
        segmentation_model,
        output_file=output_filename,
        suffix="segm",
        force=True,
    )
