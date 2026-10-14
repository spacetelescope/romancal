"""Tests for HEALPix sky offset solver and parquet offset I/O."""

from __future__ import annotations

import numpy as np
import pytest
from roman_datamodels.datamodels import ImageModel

from romancal.datamodels import ModelLibrary
from romancal.skymatch.healpix_match import (
    match_healpix_skyvals,
    pipeline_level_from_epsilon,
)
from romancal.skymatch.offset_io import (
    load_offset_table,
    load_skyval_parquet,
    lookup_offset,
    normalize_exposure_key,
    resolve_offsets_path_from_asn,
    write_offset_table,
    write_skyval_parquet,
)
from romancal.skymatch.skymatch_step import SkyMatchStep


def test_normalize_exposure_key():
    """Purpose: cal/segm/path variants collapse to one offset-table key."""
    assert normalize_exposure_key("a/b/c_cal.asdf") == "c"
    assert normalize_exposure_key("c_segm.asdf") == "c"
    assert normalize_exposure_key("c") == "c"


def test_pipeline_level_sign():
    """Purpose: level is minus epsilon so downstream subtraction matches additive fix."""
    assert pipeline_level_from_epsilon(0.01) == pytest.approx(-0.01)


def test_match_healpix_two_exposures_constant_offset():
    """
    Purpose: two exposures with shared healpix cells and a constant sky delta
    recover that delta (up to overall recenter).
    """
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 1.25, dtype=np.float64)
    samples = [
        {"pixel": pix, "value": sky0, "unc": unc},
        {"pixel": pix, "value": sky1, "unc": unc},
    ]
    eps, names = match_healpix_skyvals(samples, names=["a", "b"])
    assert names == ["a", "b"]
    # After recenter, relative difference should be ~0.25
    assert (eps[0] - eps[1]) == pytest.approx(0.25, abs=1e-6)
    # Corrected values equal
    assert np.median(sky0 + eps[0]) == pytest.approx(np.median(sky1 + eps[1]), abs=1e-6)


def test_parquet_offset_roundtrip(tmp_path):
    """Purpose: offsets parquet stores level=-epsilon and loads as pipeline levels."""
    path = write_offset_table(
        tmp_path / "off.parquet",
        ["exp_a_cal.asdf", "exp_b"],
        [0.01, -0.02],
    )
    table = load_offset_table(path)
    assert lookup_offset("exp_a_cal.asdf", table) == pytest.approx(-0.01)
    assert lookup_offset("path/exp_b_segm.asdf", table) == pytest.approx(0.02)


def test_skyval_parquet_roundtrip(tmp_path):
    """Purpose: skyval parquet + cov sidecar round-trip pixel/value/unc."""
    pix = np.array([1, 2, 3], dtype=np.int64)
    val = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    unc = np.array([0.01, 0.01, 0.01], dtype=np.float32)
    cov = np.array([10, 11], dtype=np.int64)
    out = write_skyval_parquet(
        tmp_path / "x.parquet",
        pixel=pix,
        value=val,
        unc=unc,
        coverage=cov,
        filename="x",
    )
    sample = load_skyval_parquet(out)
    np.testing.assert_array_equal(sample["pixel"], pix)
    np.testing.assert_allclose(sample["value"], val)
    np.testing.assert_allclose(sample["unc"], unc)
    np.testing.assert_array_equal(sample["coverage"], cov)


def test_resolve_offsets_path_from_asn(tmp_path):
    """Purpose: product-level ASN field resolves relative to ASN directory."""
    off = tmp_path / "offsets.parquet"
    write_offset_table(off, ["a"], [0.0])
    asn_path = tmp_path / "mos_asn.json"
    asn = {
        "products": [
            {
                "name": "prod",
                "skymatch_offsets_file": "offsets.parquet",
                "members": [],
            }
        ]
    }
    path = resolve_offsets_path_from_asn(asn, asn_path=asn_path)
    assert path == off.resolve()


def _mk_image(filename, level=0.0, shape=(20, 20)):
    model = ImageModel.create_fake_data(shape=shape)
    model.meta.filename = filename
    model.data = np.full(shape, level, dtype=np.float32)
    model.meta["background"] = {"level": None, "subtracted": False, "method": None}
    model.meta.cal_step = {}
    for step_name in model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        model.meta.cal_step[step_name] = "INCOMPLETE"
    model.meta.cal_logs = []
    return model


def test_skymatch_step_offset_file(tmp_path):
    """
    Purpose: SkyMatchStep offset_file applies parquet levels without stcal match.
    """
    off = write_offset_table(
        tmp_path / "offsets.parquet",
        ["img_a", "img_b"],
        [0.05, -0.05],  # epsilon
    )
    # levels stored as -epsilon => -0.05 and +0.05
    lib = ModelLibrary(
        [
            _mk_image("img_a_cal.asdf", 10.0),
            _mk_image("img_b_cal.asdf", 10.0),
        ]
    )
    result = SkyMatchStep.call(lib, offset_file=str(off), subtract=False)
    with result:
        models = list(result)
        assert models[0].meta.background.method == "file"
        assert models[0].meta.background.level == pytest.approx(-0.05)
        assert models[1].meta.background.level == pytest.approx(0.05)
        assert models[0].meta.cal_step.skymatch == "COMPLETE"
        for m in models:
            result.shelve(m, modify=False)


def test_skymatch_step_offset_file_missing_member(tmp_path):
    """Purpose: missing offset row raises a clear error."""
    off = write_offset_table(tmp_path / "offsets.parquet", ["img_a"], [0.0])
    lib = ModelLibrary(
        [
            _mk_image("img_a_cal.asdf"),
            _mk_image("img_b_cal.asdf"),
        ]
    )
    with pytest.raises(ValueError, match="No offset entry"):
        SkyMatchStep.call(lib, offset_file=str(off))
