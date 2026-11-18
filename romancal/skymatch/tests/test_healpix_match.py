"""Tests for HEALPix sky level solver and parquet offset I/O."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from roman_datamodels.datamodels import ImageModel

from romancal.datamodels import ModelLibrary
from romancal.skymatch.healpix_match import match_healpix_skyvals
from romancal.skymatch.offset_io import (
    load_offset_table,
    load_skyval_parquet,
    lookup_offset,
    normalize_exposure_key,
    offsets_path_from_asn,
    write_offset_table,
    write_skyval_parquet,
)
from romancal.skymatch.skymatch_step import SkyMatchStep


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("a/b/c_cal.asdf", "c"),
        ("c_segm.asdf", "c"),
        ("c", "c"),
        ("test_cal.fits", "test"),
        ("path/to/test_segm.asdf", "test"),
    ],
)
def test_normalize_exposure_key(name, expected):
    """Purpose: cal/segm/path variants collapse to one offset-table key."""
    assert normalize_exposure_key(name) == expected


def test_match_healpix_two_exposures_constant_offset():
    """
    Purpose: two exposures with shared healpix cells and a constant sky delta
    recover that delta as subtractive levels (up to overall recenter).
    """
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 1.25, dtype=np.float64)
    samples = [
        {"pixel": pix, "value": sky0, "unc": unc},
        {"pixel": pix, "value": sky1, "unc": unc},
    ]
    levels, names = match_healpix_skyvals(samples, names=["a", "b"])
    assert names == ["a", "b"]
    # corrected_value = value - level; relative delta 0.25
    assert (levels[1] - levels[0]) == pytest.approx(0.25, abs=1e-5)
    assert np.median(sky0 - levels[0]) == pytest.approx(
        np.median(sky1 - levels[1]), abs=1e-5
    )


def test_match_healpix_recenter_changes_solution():
    """Purpose: recenter=True vs False can change absolute levels (3-exp case)."""
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    # Asymmetric sky values so biweight center differs from sum(eps)=0 alone
    # when prior is disabled.
    samples = [
        {"pixel": pix, "value": np.full(pix.shape, 1.0), "unc": unc},
        {"pixel": pix, "value": np.full(pix.shape, 1.1), "unc": unc},
        {"pixel": pix, "value": np.full(pix.shape, 3.0), "unc": unc},
    ]
    lev_on, _ = match_healpix_skyvals(
        samples, names=["a", "b", "c"], recenter=True, prior_sigma_dn=0.0
    )
    lev_off, _ = match_healpix_skyvals(
        samples, names=["a", "b", "c"], recenter=False, prior_sigma_dn=0.0
    )
    assert not np.allclose(lev_on, lev_off, atol=1e-8)
    # pairwise relative differences are gauge-invariant
    assert (lev_on[1] - lev_on[0]) == pytest.approx(lev_off[1] - lev_off[0], abs=1e-6)


def test_match_healpix_min_overlap_cells():
    """Purpose: pairs below min_overlap_cells are ignored."""
    pix = np.arange(100, dtype=np.int64)  # 100 shared cells
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 1.25, dtype=np.float64)
    samples = [
        {"pixel": pix, "value": sky0, "unc": unc},
        {"pixel": pix, "value": sky1, "unc": unc},
    ]
    # threshold above overlap -> no pairs; prior pulls levels ~ 0
    lev_hi, _ = match_healpix_skyvals(
        samples, names=["a", "b"], min_overlap_cells=101, recenter=False
    )
    assert np.allclose(lev_hi, 0.0, atol=1e-6)
    lev_ok, _ = match_healpix_skyvals(
        samples, names=["a", "b"], min_overlap_cells=100, recenter=False
    )
    assert abs(lev_ok[1] - lev_ok[0]) == pytest.approx(0.25, abs=1e-4)


def test_match_healpix_photom_scale():
    """Purpose: conversion_megajanskys scales values before solve."""
    pix = np.arange(100, 200, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky0 = np.full(pix.shape, 1.0, dtype=np.float64)
    sky1 = np.full(pix.shape, 1.25, dtype=np.float64)
    conv = 2.0
    samples = [
        {
            "pixel": pix,
            "value": sky0,
            "unc": unc,
            "conversion_megajanskys": conv,
        },
        {
            "pixel": pix,
            "value": sky1,
            "unc": unc,
            "conversion_megajanskys": conv,
        },
    ]
    levels, _ = match_healpix_skyvals(samples, names=["a", "b"])
    # delta in DN/s is 0.25 -> 0.5 after *2
    assert (levels[1] - levels[0]) == pytest.approx(0.5, abs=1e-4)


def test_match_healpix_no_overlap_prior_zero():
    """Purpose: no overlaps returns near-zero levels via prior."""
    pix0 = np.array([0, 1], dtype=np.int64)
    pix1 = np.array([100, 101], dtype=np.int64)
    unc = np.full(2, 0.01, dtype=np.float64)
    samples = [
        {"pixel": pix0, "value": np.ones(2), "unc": unc},
        {"pixel": pix1, "value": np.full(2, 2.0), "unc": unc},
    ]
    levels, names = match_healpix_skyvals(samples, names=["a", "b"])
    assert names == ["a", "b"]
    assert np.allclose(levels, 0.0, atol=1e-6)


def test_match_healpix_all_bad_samples_zero_levels():
    """Purpose: all-bad samples yield zero levels (no usable pairs)."""
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
    levels, names = match_healpix_skyvals([bad, good], names=["bad", "good"])
    assert names == ["bad", "good"]
    assert np.allclose(levels, 0.0, atol=1e-5)


def test_parquet_offset_roundtrip(tmp_path):
    """Purpose: offsets parquet stores level and loads as pipeline levels."""
    path = write_offset_table(
        tmp_path / "off.parquet",
        ["exp_a_cal.asdf", "exp_b"],
        [-0.01, 0.02],
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
        conversion_megajanskys=1.5,
    )
    sample = load_skyval_parquet(out)
    np.testing.assert_array_equal(sample["pixel"], pix)
    np.testing.assert_allclose(sample["value"], val)
    np.testing.assert_allclose(sample["unc"], unc)
    np.testing.assert_array_equal(sample["coverage"], cov)
    assert sample["conversion_megajanskys"] == pytest.approx(1.5)


def test_offsets_path_from_asn_top_level():
    """Purpose: top-level ASN field is returned as-is (path join via library)."""
    asn = {"skymatch_offsets_file": "offsets.parquet"}
    assert offsets_path_from_asn(asn) == "offsets.parquet"


def test_offsets_path_from_asn_ignores_product_level():
    """Purpose: product-level offsets field is ignored (top-level only)."""
    asn = {
        "products": [
            {"name": "prod", "skymatch_offsets_file": "offsets.parquet", "members": []}
        ]
    }
    assert offsets_path_from_asn(asn) is None


def _mk_image(filename, level=0.0, shape=(20, 20), subtracted=False):
    model = ImageModel.create_fake_data(shape=shape)
    model.meta.filename = filename
    model.data = np.full(shape, level, dtype=np.float32)
    model.meta["background"] = {
        "level": None,
        "subtracted": subtracted,
        "method": "None",
    }
    model.meta.cal_step = {}
    for step_name in model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        model.meta.cal_step[step_name] = "INCOMPLETE"
    model.meta.cal_logs = []
    return model


def test_library_resolves_relative_offset_path(tmp_path):
    """Purpose: ModelLibrary._member_to_filename joins relative offset paths."""
    off = tmp_path / "offsets.parquet"
    write_offset_table(off, ["a"], [0.0])
    img = _mk_image("img_a_cal.asdf")
    img.save(tmp_path / "img_a_cal.asdf")
    asn = {
        "asn_type": "image",
        "asn_rule": "candidate",
        "version_id": "1",
        "code_version": "1",
        "target": "t",
        "program": "p",
        "asn_id": "a",
        "asn_pool": "none",
        "skymatch_offsets_file": "offsets.parquet",
        "products": [
            {
                "name": "prod",
                "members": [
                    {"expname": "img_a_cal.asdf", "exptype": "science"},
                ],
            }
        ],
    }
    asn_path = tmp_path / "asn.json"
    asn_path.write_text(json.dumps(asn))
    lib = ModelLibrary(str(asn_path))
    resolved = lib._member_to_filename({"expname": "offsets.parquet"})
    assert Path(resolved).resolve() == off.resolve()


def test_skymatch_step_offset_file(tmp_path):
    """
    Purpose: skymethod=offset_file applies parquet levels without stcal match.
    """
    off = write_offset_table(
        tmp_path / "offsets.parquet",
        ["img_a", "img_b"],
        [-0.05, 0.05],
    )
    lib = ModelLibrary(
        [
            _mk_image("img_a_cal.asdf", 10.0),
            _mk_image("img_b_cal.asdf", 10.0),
        ]
    )
    result = SkyMatchStep.call(
        lib, skymethod="offset_file", offset_file=str(off), subtract=False
    )
    with result:
        models = list(result)
        assert models[0].meta.background.method == "offset_file"
        assert models[0].meta.background.level == pytest.approx(-0.05)
        assert models[1].meta.background.level == pytest.approx(0.05)
        assert models[0].meta.background["offset_file"] == str(off)
        assert models[0].meta.cal_step.skymatch == "COMPLETE"
        for m in models:
            result.shelve(m, modify=False)


def test_skymatch_step_offset_file_ignored_unless_method(tmp_path):
    """Purpose: offset_file is ignored when skymethod is not offset_file."""
    off = write_offset_table(tmp_path / "offsets.parquet", ["img_a"], [5.0])
    lib = ModelLibrary([_mk_image("img_a_cal.asdf", 10.0)])
    # default skymethod=match should ignore offset_file (and may fail without
    # realistic WCS/data — use a library that still completes local match)
    # Just ensure we don't apply level=5 from the file when method != offset_file.
    # Running full stcal match on fake data is ok for method recording.
    result = SkyMatchStep.call(lib, offset_file=str(off), skymethod="local")
    with result:
        models = list(result)
        assert models[0].meta.background.method != "offset_file"
        assert models[0].meta.background.get("offset_file") in (None, "")
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
        SkyMatchStep.call(lib, skymethod="offset_file", offset_file=str(off))


def test_skymatch_step_offset_file_already_subtracted(tmp_path):
    """Purpose: offset_file path errors if background already subtracted."""
    off = write_offset_table(tmp_path / "offsets.parquet", ["img_a"], [0.0])
    lib = ModelLibrary([_mk_image("img_a_cal.asdf", subtracted=True)])
    with pytest.raises(ValueError, match="already subtracted"):
        SkyMatchStep.call(lib, skymethod="offset_file", offset_file=str(off))


def test_skymatch_step_subtract_changes_data(tmp_path):
    """Purpose: subtract=True subtracts pipeline level from image data."""
    off = write_offset_table(tmp_path / "offsets.parquet", ["img_a"], [-0.05])
    lib = ModelLibrary([_mk_image("img_a_cal.asdf", 10.0)])
    result = SkyMatchStep.call(
        lib, skymethod="offset_file", offset_file=str(off), subtract=True
    )
    with result:
        models = list(result)
        assert models[0].meta.background.level == pytest.approx(-0.05)
        assert np.mean(models[0].data) == pytest.approx(10.0 - (-0.05), abs=1e-5)
        for m in models:
            result.shelve(m, modify=False)


def test_skymatch_step_offset_file_missing_path(tmp_path):
    """Purpose: missing offset_file path raises FileNotFoundError."""
    lib = ModelLibrary([_mk_image("img_a_cal.asdf")])
    with pytest.raises(FileNotFoundError, match="offset_file not found"):
        SkyMatchStep.call(
            lib,
            skymethod="offset_file",
            offset_file=str(tmp_path / "nope.parquet"),
        )


def test_load_skyval_missing_columns(tmp_path):
    """Purpose: load_skyval_parquet raises for missing required columns."""
    path = tmp_path / "bad.parquet"
    Table({"healpix17": [1]}).write(path, format="parquet", overwrite=True)
    with pytest.raises(ValueError, match="missing required column"):
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
    with pytest.raises(ValueError, match="missing level column"):
        load_offset_table(path)


def test_write_offset_length_mismatch(tmp_path):
    """Purpose: write_offset_table rejects mismatched filenames/levels."""
    with pytest.raises(ValueError, match="same length"):
        write_offset_table(tmp_path / "off.parquet", ["a"], [0.1, 0.2])


def test_lookup_not_found():
    """Purpose: lookup_offset raises KeyError when not found."""
    table = {"exp_a": -0.01}
    with pytest.raises(KeyError, match="No skymatch offset"):
        lookup_offset("exp_b.asdf", table)


def test_match_healpix_empty_samples():
    """Purpose: empty samples list returns empty levels and names."""
    levels, names = match_healpix_skyvals([])
    assert len(levels) == 0
    assert names == []


def test_match_healpix_names_mismatch():
    """Purpose: names length mismatch raises ValueError."""
    pix = np.arange(100, dtype=np.int64)
    unc = np.full(pix.shape, 0.01, dtype=np.float64)
    sky = np.full(pix.shape, 1.0, dtype=np.float64)
    samples = [{"pixel": pix, "value": sky, "unc": unc}]
    with pytest.raises(ValueError, match="names must match"):
        match_healpix_skyvals(samples, names=["a", "b"])
