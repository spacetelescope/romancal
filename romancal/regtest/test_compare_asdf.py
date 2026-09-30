"""Unit tests for compare_asdf; these do not require --bigdata."""

import asdf
import numpy as np
import pytest
from astropy.table import Table

from romancal.regtest.regtestdata import compare_asdf


def _table(names=("a", "b", "c"), ee=(0.5, 0.8)):
    flux = np.arange(len(names), dtype="f8")
    table = Table({"type": list(names), "flux": flux})
    table.meta["ee_fractions"] = {"f158": np.array(ee)}
    return table


def _write(path, table):
    tree = {"roman": {"data": np.zeros(3), "sources": table}}
    asdf.AsdfFile(tree).write_to(path)
    return path


@pytest.fixture
def truth(tmp_path):
    return _write(tmp_path / "truth.asdf", _table())


def test_identical_tables(tmp_path, truth):
    result = _write(tmp_path / "result.asdf", _table())
    assert compare_asdf(result, truth).identical


def test_string_column_differs(tmp_path, truth):
    result = _write(tmp_path / "result.asdf", _table(names=("a", "x", "c")))
    diff = compare_asdf(result, truth)
    assert not diff.identical
    assert "type" in str(diff.report())


def test_table_meta_array_differs(tmp_path, truth):
    result = _write(tmp_path / "result.asdf", _table(ee=(0.5, 0.9)))
    diff = compare_asdf(result, truth)
    assert not diff.identical
    assert "metas_differ" in str(diff.report())


def test_ignored_table_not_compared(tmp_path, truth):
    # the operators must not run on ignored paths at all, not just have
    # their results dropped
    result = _write(tmp_path / "result.asdf", _table(names=("x", "y", "z")))
    diff = compare_asdf(result, truth, ignore=["roman.sources"])
    assert diff.identical
