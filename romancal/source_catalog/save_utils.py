import logging

import numpy as np
from photutils.segmentation import SegmentationImage
from roman_datamodels.datamodels import (
    ForcedImageSourceCatalogModel,
    ForcedMosaicSourceCatalogModel,
    ImageModel,
    ImageSourceCatalogModel,
    MosaicSegmentationMapModel,
    MosaicSourceCatalogModel,
    MultibandSegmentationMapModel,
    SegmentationMapModel,
)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_segment_image(
    self, segment_img, source_catalog_model, output_filename, save_debug_info=False
):
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
        log.warning("No segmentation image to save.")
        return

    # Define the segmentation model based on the source catalog model
    # type
    if isinstance(
        source_catalog_model, (ImageSourceCatalogModel | ForcedImageSourceCatalogModel)
    ):
        segm_model = SegmentationMapModel
    elif isinstance(
        source_catalog_model,
        (MosaicSourceCatalogModel | ForcedMosaicSourceCatalogModel),
    ):
        segm_model = MosaicSegmentationMapModel
    else:
        segm_model = MultibandSegmentationMapModel
    segmentation_model = segm_model.create_minimal({"meta": source_catalog_model.meta})

    # carry over image_metas if it exists (since it's not required in the schemas)
    if image_metas := source_catalog_model.meta.get("image_metas"):
        segmentation_model.meta.image_metas = image_metas

    # Set the data and detection image
    segmentation_model.data = segment_img.data.astype(np.uint32)
    if hasattr(segment_img, "detection_image"):
        segmentation_model["detection_image"] = segment_img.detection_image

    # Source injection data
    if hasattr(segment_img, "injected_sources"):
        segmentation_model["injected_sources"] = segment_img.injected_sources
    if hasattr(segment_img, "recovered_sources"):
        segmentation_model["recovered_sources"] = segment_img.recovered_sources

    # Write data for tests
    if save_debug_info:
        if hasattr(segment_img, "si_segment_img"):
            segmentation_model["si_data"] = segment_img.si_segment_img.data.astype(
                np.uint32
            )
        if hasattr(segment_img, "si_detection_image"):
            segmentation_model["si_detection_image"] = segment_img.si_detection_image

    # Save the segmentation image to the output file
    self.output_ext = "asdf"
    self.save_model(
        segmentation_model,
        output_file=output_filename,
        suffix="segm",
        force=True,
    )


def save_all_results(
    self, segment_img, cat_model, input_model=None, save_debug_info=False
):
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
    self : `romancal.source_catalog.SourceCatalogStep` or `romancal.multiband_catalog.MultiBandCatalogStep`
        The step instance that is calling this function.

    segment_img : `photutils.segmentation.SegmentationImage`
        The segmentation image created by the source catalog step.

    cat_model : `datamodels.ImageSourceCatalogModel`, `datamodels.MosaicSourceCatalogModel`, `datamodels.ForcedImageSourceCatalogModel`, or `datamodels.ForcedMosaicSourceCatalogModel`
        The source catalog model created by the step

    input_model : `None`, `ImageModel`, or `MosaicModel`, optional
        The input model to the source catalog step. This
        is required only for the SourceCatalogStep when
        ``self.return_updated_model`` is `True`.

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
        self.output_file if self.output_file is not None else cat_model.meta.filename
    )

    # always save the segmentation image
    save_segment_image(self, segment_img, cat_model, output_filename, save_debug_info)

    # Update the source catalog filename metadata
    self.output_ext = "parquet"
    output_catalog_name = self.make_output_path(basepath=output_filename, suffix="cat")
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
    if return_updated_model:
        # define the catalog filename; self.save_model will
        # determine whether to use a fully qualified path

        # set the suffix to something else to prevent the step from
        # overwriting the source catalog file with a datamodel
        self.suffix = "sourcecatalog"

        # update the input datamodel with the tweakreg catalog name
        if isinstance(input_model, ImageModel):
            input_model.meta.source_catalog = {
                "tweakreg_catalog_name": output_catalog_name
            }
            input_model.meta.cal_step.source_catalog = "COMPLETE"

        result = input_model
    else:
        self.output_ext = "parquet"
        result = cat_model

    if self.output_file is None:
        self.output_file = cat_model.meta.filename

    # validate the result to flush out any lazy-loaded contents
    result.validate()
    return result


def save_empty_results(self, segment_shape, cat_model, input_model=None, msg=""):
    """
    Return and save empty results for the source catalog or multiband
    catalog step.

    Parameters
    ----------
    self : `romancal.source_catalog.SourceCatalogStep` or `romancal.multiband_catalog.MultiBandCatalogStep`
        The step instance that is calling this function.

    segment_shape : tuple
        The shape of the segmentation image to create.

    cat_model : `datamodels.ImageSourceCatalogModel`, `datamodels.MosaicSourceCatalogModel`, `datamodels.ForcedImageSourceCatalogModel`, or `datamodels.ForcedMosaicSourceCatalogModel`
        The source catalog model created by the step.

    input_model : `None`, `ImageModel`, or `MosaicModel`, optional
        The input model to the source catalog step. This
        is required only for the SourceCatalogStep when
        ``self.return_updated_model`` is `True`.

    msg : str, optional
        A message to log if no results were found. Default is an empty
        string.

    Returns
    -------
    result : `ImageModel`, `MosaicModel`, or `SourceCatalog`
        The source catalog model or the input model with updated
        metadata.
    """
    if msg:
        log.error(msg)

    segment_img = np.zeros(segment_shape, dtype=np.uint32)
    cat_model.source_catalog = cat_model.create_empty_catalog()
    return save_all_results(self, segment_img, cat_model, input_model=input_model)
