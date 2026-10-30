"""Tests for HEALPix sky offset solver and parquet offset I/O."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.table import Table
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


def test_match_healpix_empty_samples():
    """Purpose: empty samples list returns empty epsilon and names."""
    eps, names = match_healpix_skyvals([])
    assert len(eps) == 0
    assert names == []


def test_match_healpix_default_names():
    """Purpose: default names are exp000, exp001, etc."""
    pix = np.arange(100, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky = np.full(pix.shape, 1.0, dtype=np.float64)
    samples = [
        {"pixel": pix, "value": sky, "unc": unc},
        {"pixel": pix, "value": sky, "unc": unc},
    ]
    eps, names = match_healpix_skyvals(samples)
    assert names == ["exp000", "exp001"]


def test_match_healpix_names_mismatch():
    """Purpose: names length mismatch raises ValueError."""
    pix = np.arange(100, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky = np.full(pix.shape, 1.0, dtype=np.float64)
    samples = [{"pixel": pix, "value": sky, "unc": unc}]
    with pytest.raises(ValueError, match="names must match"):
        match_healpix_skyvals(samples, names=["a", "b"])


def test_match_healpix_length_mismatch():
    """Purpose: pixel/value/unc length mismatch raises ValueError."""
    pix = np.arange(100, dtype=np.int64)
    unc = np.full(50, 0.01, dtype=np.float64)
    sky = np.full(100, 1.0, dtype=np.float64)
    samples = [{"pixel": pix, "value": sky, "unc": unc}]
    with pytest.raises(ValueError, match="pixel/value/unc length mismatch"):
        match_healpix_skyvals(samples)


def test_match_healpix_bad_pixels():
    """Purpose: non-finite and non-positive uncertainties are filtered."""
    pix = np.array([0, 1, 2, 3, 4], dtype=np.int64)
    unc = np.array([0.01, np.nan, -0.01, 0.0, 0.02], dtype=np.float64)
    sky = np.full(5, 1.0, dtype=np.float64)
    samples = [{"pixel": pix, "value": sky, "unc": unc}]
    eps, names = match_healpix_skyvals(samples, names=["a"])
    assert len(names) == 1


def test_match_healpix_no_overlap():
    """Purpose: no overlapping HEALPix cells returns zeros."""
    pix0 = np.array([0, 1], dtype=np.int64)
    pix1 = np.array([100, 101], dtype=np.int64)
    unc = np.full(2, 0.01, dtype=np.float64)
    sky0 = np.full(2, 1.0, dtype=np.float64)
    sky1 = np.full(2, 2.0, dtype=np.float64)
    samples = [
        {"pixel": pix0, "value": sky0, "unc": unc},
        {"pixel": pix1, "value": sky1, "unc": unc},
    ]
    eps, names = match_healpix_skyvals(samples, names=["a", "b"])
    assert len(eps) == 2
    assert all(eps == 0)


def test_match_healpix_coverage_no_overlap():
    """Purpose: coarse coverage skips pairs with no overlap."""
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 2.0, dtype=np.float64)
    cov0 = np.array([0, 1], dtype=np.int64)
    cov1 = np.array([2, 3], dtype=np.int64)
    samples = [
        {"pixel": pix, "value": sky0, "unc": unc, "coverage": cov0},
        {"pixel": pix, "value": sky1, "unc": unc, "coverage": cov1},
    ]
    eps, names = match_healpix_skyvals(samples, names=["a", "b"])
    assert len(eps) == 2


def test_match_healpix_recenter_false():
    """Purpose: recenter=False does not apply biweight centering."""
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 1.25, dtype=np.float64)
    samples = [
        {"pixel": pix, "value": sky0, "unc": unc},
        {"pixel": pix, "value": sky1, "unc": unc},
    ]
    eps, names = match_healpix_skyvals(samples, names=["a", "b"], recenter=False)
    assert len(eps) == 2


def test_match_healpix_min_overlap():
    """Purpose: min_overlap_cells filters pairs below threshold."""
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 1.25, dtype=np.float64)
    samples = [
        {"pixel": pix, "value": sky0, "unc": unc},
        {"pixel": pix, "value": sky1, "unc": unc},
    ]
    eps, names = match_healpix_skyvals(
        samples, names=["a", "b"], min_overlap_cells=1000
    )
    assert len(eps) == 2


def test_pipeline_level_array():
    """Purpose: pipeline_level_from_epsilon works with arrays."""
    vals = pipeline_level_from_epsilon(np.array([0.01, 0.02]))
    np.testing.assert_allclose(vals, [-0.01, -0.02])


def test_offset_as_pipeline_level(tmp_path):
    """Purpose: write_offset_table with as_pipeline_level=True."""
    path = write_offset_table(
        tmp_path / "off.parquet", ["a"], [0.01], as_pipeline_level=True
    )
    table = load_offset_table(path)
    # When as_pipeline_level=True: level=0.01, epsilon=-0.01
    assert lookup_offset("a", table) == pytest.approx(0.01)


def test_offset_top_level_asn():
    """Purpose: resolve_offsets_path_from_asn finds ASN-level field."""
    asn = {"skymatch_offsets_file": "offsets.parquet"}
    path = resolve_offsets_path_from_asn(asn)
    assert path is not None
    assert path.name == "offsets.parquet"


def test_offset_no_candidates():
    """Purpose: resolve_offsets_path_from_asn returns None when no field."""
    asn = {"products": []}
    path = resolve_offsets_path_from_asn(asn)
    assert path is None


def test_lookup_raw_fallback():
    """Purpose: lookup_offset falls back to raw stem when key not found."""
    table = {"exp_a": -0.01}
    result = lookup_offset("path/exp_a.asdf", table)
    assert result == pytest.approx(-0.01)


def test_lookup_not_found():
    """Purpose: lookup_offset raises KeyError when not found."""
    table = {"exp_a": -0.01}
    with pytest.raises(KeyError, match="No skymatch offset"):
        lookup_offset("exp_b.asdf", table)


def test_skyval_no_coverage(tmp_path):
    """Purpose: skyval parquet works without coverage sidecar."""
    pix = np.array([1, 2, 3], dtype=np.int64)
    val = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    unc = np.array([0.01, 0.01, 0.01], dtype=np.float32)
    out = write_skyval_parquet(tmp_path / "x.parquet", pixel=pix, value=val, unc=unc)
    sample = load_skyval_parquet(out)
    assert "coverage" not in sample


def test_normalize_extra():
    """Purpose: additional normalize_exposure_key variants."""
    assert normalize_exposure_key("test_cal.fits") == "test"
    assert normalize_exposure_key("test_i2d.asdf") == "test"
    assert normalize_exposure_key("path/to/test_segm.asdf") == "test"


def test_write_offset_suffix(tmp_path):
    """Purpose: write_offset_table adds .parquet suffix when missing."""
    path = write_offset_table(tmp_path / "off", ["a"], [0.0])
    assert path.suffix == ".parquet"


def test_load_offset_filename_column(tmp_path):
    """Purpose: load_offset_table handles filename column."""
    table = Table({"filename": ["exp1", "exp2"], "level": [0.1, -0.1]})
    path = tmp_path / "test.parquet"
    table.write(path, format="parquet", overwrite=True)
    result = load_offset_table(path)
    assert "exp1" in result


def test_load_offset_file_column(tmp_path):
    """Purpose: load_offset_table handles 'file' column."""
    table = Table({"file": ["exp1"], "level": [0.1]})
    path = tmp_path / "test.parquet"
    table.write(path, format="parquet", overwrite=True)
    result = load_offset_table(path)
    assert "exp1" in result


def test_load_offset_alt_columns(tmp_path):
    """Purpose: load_offset_table handles epsilon/offset/offsets columns."""
    for col in ["epsilon", "offset", "offsets"]:
        table = Table({"filename": ["exp1"], col: [0.05]})
        path = tmp_path / f"test_{col}.parquet"
        table.write(path, format="parquet", overwrite=True)
        result = load_offset_table(path)
        assert "exp1" in result


def test_skymatch_step_subtract(tmp_path):
    """Purpose: SkyMatchStep with subtract=True modifies data."""
    off = write_offset_table(tmp_path / "offsets.parquet", ["img_a"], [0.05])
    lib = ModelLibrary(
        [
            _mk_image("img_a_cal.asdf", 10.0),
        ]
    )
    result = SkyMatchStep.call(lib, offset_file=str(off), subtract=True)
    with result:
        models = list(result)
        assert models[0].meta.background.method == "file"
        for m in models:
            result.shelve(m, modify=False)


def test_lookup_raw_stem_fallback():
    """Purpose: lookup uses Path.stem when normalized key is absent."""
    # normalize("x_cal.asdf") -> "x"; stem remains "x_cal"
    table = {"x_cal": -0.01}
    assert lookup_offset("x_cal.asdf", table) == pytest.approx(-0.01)


def test_write_skyval_adds_suffix(tmp_path):
    """Purpose: write_skyval_parquet appends .parquet when missing."""
    out = write_skyval_parquet(
        tmp_path / "x",
        pixel=np.array([1], dtype=np.int64),
        value=np.array([0.1], dtype=np.float32),
        unc=np.array([0.01], dtype=np.float32),
    )
    assert out.suffix == ".parquet"
    assert out.exists()


def test_load_skyval_alt_columns(tmp_path):
    """Purpose: load_skyval_parquet accepts pixel/value/unc column names."""
    path = tmp_path / "alt.parquet"
    Table(
        {
            "pixel": np.array([1, 2], dtype=np.int64),
            "value": np.array([0.1, 0.2], dtype=np.float64),
            "unc": np.array([0.01, 0.01], dtype=np.float64),
        }
    ).write(path, format="parquet", overwrite=True)
    sample = load_skyval_parquet(path)
    np.testing.assert_array_equal(sample["pixel"], [1, 2])
    np.testing.assert_allclose(sample["value"], [0.1, 0.2])
    np.testing.assert_allclose(sample["unc"], [0.01, 0.01])


def test_load_skyval_missing_columns(tmp_path):
    """Purpose: load_skyval_parquet raises for missing required columns."""
    path = tmp_path / "bad.parquet"
    Table({"healpix17": [1]}).write(path, format="parquet", overwrite=True)
    with pytest.raises(ValueError, match="missing data/value"):
        load_skyval_parquet(path)

    Table({"data": [0.1]}).write(path, format="parquet", overwrite=True)
    with pytest.raises(ValueError, match="missing healpix17/pixel"):
        load_skyval_parquet(path)

    Table({"healpix17": [1], "data": [0.1]}).write(
        path, format="parquet", overwrite=True
    )
    with pytest.raises(ValueError, match="missing err/unc"):
        load_skyval_parquet(path)


def test_load_offset_missing_columns(tmp_path):
    """Purpose: load_offset_table raises when name or level columns are absent."""
    path = tmp_path / "bad.parquet"
    Table({"level": [0.1]}).write(path, format="parquet", overwrite=True)
    with pytest.raises(ValueError, match="missing filename column"):
        load_offset_table(path)

    Table({"filename": ["a"], "foo": [1.0]}).write(
        path, format="parquet", overwrite=True
    )
    with pytest.raises(ValueError, match="missing level/epsilon/offset"):
        load_offset_table(path)


def test_write_offset_length_mismatch(tmp_path):
    """Purpose: write_offset_table rejects mismatched filenames/offsets."""
    with pytest.raises(ValueError, match="same length"):
        write_offset_table(tmp_path / "off.parquet", ["a"], [0.1, 0.2])


def test_skymatch_step_offset_file_missing_path(tmp_path):
    """Purpose: missing offset_file path raises FileNotFoundError."""
    lib = ModelLibrary([_mk_image("img_a_cal.asdf")])
    with pytest.raises(FileNotFoundError, match="offset_file not found"):
        SkyMatchStep.call(lib, offset_file=str(tmp_path / "nope.parquet"))


def test_skymatch_step_subtract_changes_data(tmp_path):
    """Purpose: subtract=True subtracts pipeline level from image data."""
    off = write_offset_table(tmp_path / "offsets.parquet", ["img_a"], [0.05])  # epsilon
    # level stored as -epsilon = -0.05
    lib = ModelLibrary([_mk_image("img_a_cal.asdf", 10.0)])
    result = SkyMatchStep.call(lib, offset_file=str(off), subtract=True)
    with result:
        models = list(result)
        assert models[0].meta.background.level == pytest.approx(-0.05)
        assert np.mean(models[0].data) == pytest.approx(10.0 - (-0.05), abs=1e-5)
        for m in models:
            result.shelve(m, modify=False)


def test_resolve_offsets_absolute_path(tmp_path):
    """Purpose: absolute ASN offset paths are returned unchanged."""
    off = tmp_path / "offsets.parquet"
    write_offset_table(off, ["a"], [0.0])
    asn = {"skymatch_offsets_file": str(off.resolve())}
    path = resolve_offsets_path_from_asn(asn)
    assert path == off.resolve()


def test_match_healpix_empty_after_bad_filter():
    """Purpose: exposures with all-bad samples are skipped in pair loops."""
    bad = {
        "pixel": np.array([1, 2], dtype=np.int64),
        "value": np.array([1.0, 1.0], dtype=np.float64),
        "unc": np.array([np.nan, 0.0], dtype=np.float64),
    }
    good = {
        "pixel": np.array([1, 2, 3], dtype=np.int64),
        "value": np.ones(3, dtype=np.float64),
        "unc": np.full(3, 0.01, dtype=np.float64),
    }
    # Cover both empty-a and empty-b skip branches
    eps1, names1 = match_healpix_skyvals([bad, good], names=["bad", "good"])
    eps2, names2 = match_healpix_skyvals([good, bad], names=["good", "bad"])
    assert names1 == ["bad", "good"] and names2 == ["good", "bad"]
    assert len(eps1) == 2 and len(eps2) == 2
