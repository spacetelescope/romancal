"""Tests for product cal_logs policy helpers (#2510)."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
import roman_datamodels.datamodels as rdm
from astropy.time import Time

from romancal.datamodels import ModelLibrary
from romancal.multiband_catalog._metadata import (
    blend_image_metadata,
    finalize_catalog_metadata,
)
from romancal.multiband_catalog._multiband_catalog import initialize_catalog_model
from romancal.pipeline import ExposurePipeline
from romancal.resample._meta_blender import MetaBlender
from romancal.stpipe.cal_logs import (
    copy_meta_without_cal_logs,
    extract_image_meta,
    get_cal_logs,
    set_product_cal_logs,
    slice_log_records,
    strip_inherited_cal_logs,
)


def test_copy_meta_without_cal_logs_drops_logs():
    """
    Purpose: Metadata copies used to seed new products must not inherit
    upstream cal_logs.
    """
    meta = {
        "filename": "a.asdf",
        "cal_logs": ["upstream"],
        "instrument": {"name": "WFI"},
    }
    copied = copy_meta_without_cal_logs(meta)
    assert "cal_logs" not in copied
    assert copied["filename"] == "a.asdf"
    assert copied["instrument"]["name"] == "WFI"
    # deep copy: mutating result does not affect input
    copied["instrument"]["name"] = "X"
    assert meta["instrument"]["name"] == "WFI"


def test_extract_image_meta_excludes_cal_logs_and_wcs():
    """
    Purpose: Multiband image_metas must omit cal_logs (and wcs) so L3
    production logs do not land in catalog parquet metadata.
    """
    meta = {
        "filename": "l3.asdf",
        "cal_logs": ["L3LOG_0", "L3LOG_1"],
        "wcs": object(),
        "individual_image_meta": {"basic": 1},
        "instrument": {"optical_element": "F184"},
    }
    image_meta = extract_image_meta(meta)
    assert "cal_logs" not in image_meta
    assert "wcs" not in image_meta
    assert "individual_image_meta" not in image_meta
    assert image_meta["instrument"]["optical_element"] == "F184"


def test_set_and_strip_product_cal_logs(base_image):
    """
    Purpose: set_product_cal_logs writes meta.cal_logs and strip removes it.
    """
    model = base_image()
    set_product_cal_logs(model, ["a", "b"])
    assert get_cal_logs(model) == ["a", "b"]
    strip_inherited_cal_logs(model)
    assert get_cal_logs(model) == []


def test_slice_log_records():
    """
    Purpose: Member-scoped ELP logs are taken as a slice of the pipeline
    log stream.
    """
    records = ["a", "b", "c", "d"]
    assert slice_log_records(records, 1, 3) == ["b", "c"]
    assert slice_log_records(records, 2) == ["c", "d"]
    assert slice_log_records(records, 10) == []


def test_multiband_image_metas_exclude_l3_cal_logs(tmp_path):
    """
    Purpose: Multiband catalog image_metas and parquet metadata must not
    contain L3 production cal_logs.
    """
    m1 = rdm.MosaicModel.create_fake_data(shape=(20, 20))
    m1.meta.filename = "f184.asdf"
    m1.meta.cal_logs = [f"L3A_{i}" for i in range(5)]
    m1.meta.instrument.optical_element = "F184"

    m2 = rdm.MosaicModel.create_fake_data(shape=(20, 20))
    m2.meta.filename = "f158.asdf"
    m2.meta.cal_logs = [f"L3B_{i}" for i in range(3)]
    m2.meta.instrument.optical_element = "F158"

    library = ModelLibrary([m1, m2])
    cat = initialize_catalog_model(library, m1)
    assert get_cal_logs(cat) == []

    time_means: list = []
    exposure_times: list = []
    blend_image_metadata(m1, cat, time_means, exposure_times)
    blend_image_metadata(m2, cat, time_means, exposure_times)
    finalize_catalog_metadata(cat, time_means, exposure_times)

    assert len(cat.meta.image_metas) == 2
    for image_meta in cat.meta.image_metas:
        assert "cal_logs" not in image_meta

    cat.source_catalog = cat.create_empty_catalog()
    out = tmp_path / "mb.parquet"
    cat.to_parquet(out)

    import pyarrow.parquet as pq

    meta = pq.read_schema(out).metadata or {}
    keys = [k.decode() if isinstance(k, bytes) else k for k in meta]
    cal_keys = [k for k in keys if "cal_log" in k]
    assert cal_keys == [], f"unexpected cal_logs keys in parquet: {cal_keys}"


def test_meta_blender_keeps_l2_logs_off_meta_cal_logs():
    """
    Purpose: L3 mosaics store input L2 logs on individual_image_cal_logs,
    not meta.cal_logs.
    """
    img1 = rdm.ImageModel.create_fake_data(shape=(10, 10))
    img1.meta.filename = "a.asdf"
    img1.meta.cal_logs = ["L2_A"]
    img2 = rdm.ImageModel.create_fake_data(shape=(10, 10))
    img2.meta.filename = "b.asdf"
    img2.meta.cal_logs = ["L2_B"]
    # keep exposure times comparable for blender min/max
    img2.meta.exposure.start_time = img1.meta.exposure.start_time
    img2.meta.exposure.end_time = img1.meta.exposure.end_time

    blender = MetaBlender()
    blender.blend(img1)
    blender.blend(img2)
    mosaic = blender.finalize()

    assert mosaic.meta.get("cal_logs") in (
        None,
        [],
    )
    assert mosaic["individual_image_cal_logs"] == [["L2_A"], ["L2_B"]]


def test_finalize_result_respects_member_scoped_logs(base_image):
    """
    Purpose: ExposurePipeline.finalize_result keeps member-scoped cal_logs
    instead of replacing them with the full multi-member stream.
    """
    model = base_image()
    set_product_cal_logs(model, ["member-only"])

    pipeline = ExposurePipeline()
    pipeline._member_scoped_ids = {id(model)}
    pipeline._log_records = ["full", "pipeline", "stream"]
    pipeline.finalize_result(model, [])

    assert get_cal_logs(model) == ["member-only"]


def test_elp_member_scoped_cal_logs(function_jail, tmp_path):
    """
    Purpose: Each ELP association member receives only logs emitted while
    that member was processed, not markers from sibling detectors.
    """
    paths = []
    for name in ("det1_uncal.asdf", "det2_uncal.asdf"):
        model = rdm.ScienceRawModel.create_fake_data()
        model.meta.filename = name
        model.meta.exposure.read_pattern = [[1], [2], [3], [4]]
        model.meta.exposure.ma_table_number = 5
        model.meta.exposure.start_time = Time(
            "2024-01-03T00:00:00.0", format="isot", scale="utc"
        )
        model.data = np.zeros((4, 64, 64), dtype=model.data.dtype)
        model.amp33 = np.zeros((4, 64, 128), dtype=model.amp33.dtype)
        path = tmp_path / name
        model.save(path)
        paths.append(str(path.name))

    from romancal.associations.asn_from_list import asn_from_list

    asn = asn_from_list(paths, product_name="elp_out")
    base_fn, contents = asn.dump()
    asn_path = tmp_path / base_fn
    asn_path.write_text(contents)

    pipeline = ExposurePipeline()
    pipeline.prefetch_references = False
    # Only dq_init needs to run to convert ScienceRaw -> Ramp/Image path start;
    # skip the rest for speed. Member log scoping still applies around the loop.
    for step_name in pipeline.step_defs:
        getattr(pipeline, step_name).skip = True
    pipeline.dq_init.skip = False

    # Inject distinctive log markers around each member via a thin wrapper.
    original_process_model = pipeline._process_model

    def process_model_with_marker(model):
        logging.getLogger("romancal").info("MEMBER_MARKER %s", model.meta.filename)
        return original_process_model(model)

    pipeline._process_model = process_model_with_marker

    # Avoid tweakreg work on empty products
    pipeline.tweakreg.skip = True

    # run from tmp_path so relative asn members resolve
    cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        lib, _cats, _segs = pipeline.run(str(asn_path.name))
    finally:
        os.chdir(cwd)

    with lib:
        models = list(lib)
        for m in models:
            lib.shelve(m, modify=False)

    assert len(models) == 2
    logs0 = "\n".join(get_cal_logs(models[0]))
    logs1 = "\n".join(get_cal_logs(models[1]))

    assert "MEMBER_MARKER det1_uncal.asdf" in logs0
    assert "MEMBER_MARKER det2_uncal.asdf" not in logs0
    assert "MEMBER_MARKER det2_uncal.asdf" in logs1
    assert "MEMBER_MARKER det1_uncal.asdf" not in logs1
