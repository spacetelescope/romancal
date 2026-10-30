"""Tests for roman_export_skyvals_parquet CLI helpers and main()."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from roman_datamodels.datamodels import SegmentationMapModel

from romancal.scripts.export_skyvals_parquet import (
    _expand,
    _stem_from_segm,
    main,
)
from romancal.skymatch.offset_io import load_offset_table, load_skyval_parquet


def _write_segm_with_skyvals(path: Path, *, with_cov: bool = True, n: int = 5) -> Path:
    """Create a minimal L2 segmentation product with skyvals."""
    model = SegmentationMapModel.create_fake_data(shape=(8, 8))
    dt = np.dtype(
        [
            ("healpix17", np.int64),
            ("data", np.float32),
            ("err", np.float32),
            ("covfrac", np.float32),
        ]
    )
    sky = np.zeros(n, dtype=dt)
    sky["healpix17"] = np.arange(n, dtype=np.int64)
    sky["data"] = np.linspace(1.0, 1.0 + 0.1 * (n - 1), n, dtype=np.float32)
    sky["err"] = np.full(n, 0.01, dtype=np.float32)
    sky["covfrac"] = np.ones(n, dtype=np.float32)
    model.skyvals = sky
    if with_cov:
        model.healpix11_cov = np.array([10, 11], dtype=np.int64)
    model.save(path)
    return path


def test_expand_file_glob_and_dir(tmp_path):
    """Purpose: _expand handles plain files, globs, directories, and dedupes."""
    a = tmp_path / "a_segm.asdf"
    b = tmp_path / "b_segm.asdf"
    a.touch()
    b.touch()
    (tmp_path / "c_cal.asdf").touch()

    assert _expand([str(a)]) == [a]
    globs = _expand([str(tmp_path / "*_segm.asdf")])
    assert {p.name for p in globs} == {"a_segm.asdf", "b_segm.asdf"}
    dirs = _expand([str(tmp_path)])
    assert {p.name for p in dirs} == {"a_segm.asdf", "b_segm.asdf"}
    # dedupe when same file listed multiple ways
    assert len(_expand([str(a), str(tmp_path / "a_segm.asdf")])) == 1


def test_stem_from_segm():
    """Purpose: _stem_from_segm strips trailing _segm only."""
    assert _stem_from_segm(Path("foo_segm.asdf")) == "foo"
    assert _stem_from_segm(Path("foo.asdf")) == "foo"
    assert _stem_from_segm(Path("/x/y/bar_segm.asdf")) == "bar"


def test_main_no_inputs(tmp_path):
    """Purpose: main exits when expansion finds no files."""
    with pytest.raises(SystemExit, match="No input segmentation files found"):
        main(["-o", str(tmp_path / "out"), str(tmp_path / "no_match_*.asdf")])


def test_main_missing_skyvals(tmp_path):
    """Purpose: main exits when a segm product lacks skyvals."""
    path = tmp_path / "empty_segm.asdf"
    SegmentationMapModel.create_fake_data(shape=(8, 8)).save(path)
    with pytest.raises(SystemExit, match="has no skyvals"):
        main(["-o", str(tmp_path / "out"), str(path)])


def test_main_export_and_match(tmp_path, capsys):
    """Purpose: main writes skyval parquet and optional offsets parquet."""
    seg1 = _write_segm_with_skyvals(tmp_path / "exp1_segm.asdf")
    # second exposure with offset sky values on same pixels
    seg2 = tmp_path / "exp2_segm.asdf"
    model = SegmentationMapModel.create_fake_data(shape=(8, 8))
    dt = np.dtype(
        [
            ("healpix17", np.int64),
            ("data", np.float32),
            ("err", np.float32),
            ("covfrac", np.float32),
        ]
    )
    sky = np.zeros(5, dtype=dt)
    sky["healpix17"] = np.arange(5, dtype=np.int64)
    sky["data"] = np.linspace(1.25, 1.25 + 0.1 * 4, 5, dtype=np.float32)
    sky["err"] = np.full(5, 0.01, dtype=np.float32)
    sky["covfrac"] = 1.0
    model.skyvals = sky
    model.healpix11_cov = np.array([10, 11], dtype=np.int64)
    model.save(seg2)

    outdir = tmp_path / "skyvals"
    offsets = tmp_path / "offsets.parquet"
    main(
        [
            "-o",
            str(outdir),
            "--match-offsets",
            str(offsets),
            "-v",
            str(seg1),
            str(seg2),
        ]
    )
    captured = capsys.readouterr()
    assert "Wrote 2 skyval parquet" in captured.out
    assert offsets.exists() or Path(str(offsets)).with_suffix(".parquet").exists()

    sample = load_skyval_parquet(outdir / "exp1.parquet")
    assert "pixel" in sample and "coverage" in sample
    table = load_offset_table(offsets)
    assert set(table) == {"exp1", "exp2"}


def test_main_export_without_cov_or_match(tmp_path, capsys):
    """Purpose: main works without healpix11_cov and without --match-offsets."""
    seg = _write_segm_with_skyvals(tmp_path / "only_segm.asdf", with_cov=False)
    outdir = tmp_path / "out"
    main(["-o", str(outdir), str(seg)])
    captured = capsys.readouterr()
    assert "Wrote 1 skyval parquet" in captured.out
    sample = load_skyval_parquet(outdir / "only.parquet")
    assert "coverage" not in sample
