"""Tests for product calibration log policy."""

from __future__ import annotations

import logging

import roman_datamodels.datamodels as rdm

from romancal.resample._meta_blender import MetaBlender
from romancal.stpipe import RomanStep


def test_meta_blender_does_not_store_individual_image_cal_logs():
    """
    Purpose: L3 mosaics must not carry broken/redundant individual_image_cal_logs;
    L3 meta.cal_logs already captures mosaic production messages.
    """
    img1 = rdm.ImageModel.create_fake_data(shape=(10, 10))
    img1.meta.filename = "a.asdf"
    img1.meta.cal_logs = ["L2_A"]
    img2 = rdm.ImageModel.create_fake_data(shape=(10, 10))
    img2.meta.filename = "b.asdf"
    img2.meta.cal_logs = ["L2_B"]
    img2.meta.exposure.start_time = img1.meta.exposure.start_time
    img2.meta.exposure.end_time = img1.meta.exposure.end_time

    blender = MetaBlender()
    blender.blend(img1)
    blender.blend(img2)
    mosaic = blender.finalize()

    assert "individual_image_cal_logs" not in mosaic
    # Blender itself does not set mosaic production logs
    assert mosaic.meta.get("cal_logs") in (None, [])


def test_finalize_result_sets_catalog_and_segmentation_cal_logs():
    """
    Purpose: Catalog and segmentation products receive the producing step's
    log_records (same stream on both), not parent-image logs alone.
    """
    logger_name = "test_catalog_cal_logs"
    logger = logging.getLogger(logger_name)

    class CatalogLoggingStep(RomanStep):
        def process(self, dataset):  # noqa: ARG002
            logger.info("CATALOG_STEP_MARKER")
            cat = rdm.ImageSourceCatalogModel.create_minimal()
            seg = rdm.SegmentationMapModel.create_minimal()
            # Simulate inherited parent logs that must be replaced by step logs
            cat.meta["cal_logs"] = ["PARENT_ONLY"]
            seg.meta["cal_logs"] = ["PARENT_ONLY"]
            return cat, seg

        @staticmethod
        def get_stpipe_loggers():
            return (*RomanStep.get_stpipe_loggers(), logger_name)

    cat, seg = CatalogLoggingStep().run(None)

    cat_logs = list(cat.meta.get("cal_logs") or [])
    seg_logs = list(seg.meta.get("cal_logs") or [])
    assert cat_logs == seg_logs
    assert any("CATALOG_STEP_MARKER" in line for line in cat_logs)
    assert not any("PARENT_ONLY" == line for line in cat_logs)


def test_finalize_result_sets_multiband_catalog_cal_logs():
    """
    Purpose: Multiband catalog/segmentation products get multiband-step logs.
    """
    logger_name = "test_multiband_cal_logs"
    logger = logging.getLogger(logger_name)

    class MultibandLoggingStep(RomanStep):
        def process(self, dataset):  # noqa: ARG002
            logger.info("MULTIBAND_STEP_MARKER")
            cat = rdm.MultibandSourceCatalogModel.create_minimal()
            seg = rdm.MultibandSegmentationMapModel.create_minimal()
            return cat, seg

        @staticmethod
        def get_stpipe_loggers():
            return (*RomanStep.get_stpipe_loggers(), logger_name)

    cat, seg = MultibandLoggingStep().run(None)
    assert any("MULTIBAND_STEP_MARKER" in line for line in cat.meta.cal_logs)
    assert list(cat.meta.cal_logs) == list(seg.meta.cal_logs)
