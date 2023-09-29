import asdf
import astropy.table
import numpy as np
import pytest
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.regtest.regtestdata import compare_asdf


def test_compare_asdf_identical(tmp_path):
    fn0 = tmp_path / "test0.asdf"
    fn1 = tmp_path / "test1.asdf"
    l2 = rdm.ImageModel(maker_utils.mk_level2_image(shape=(100, 100)))
    l2.save(fn0)
    l2.save(fn1)
    diff = compare_asdf(fn0, fn1)
    assert diff.identical, diff.report()


def test_compare_asdf_differ(tmp_path):
    fn0 = tmp_path / "test0.asdf"
    fn1 = tmp_path / "test1.asdf"
    l2 = rdm.ImageModel(maker_utils.mk_level2_image(shape=(100, 100)))
    l2.save(fn0)
    l2.data += 1 * l2.data.unit
    l2.save(fn1)
    diff = compare_asdf(fn0, fn1)
    assert not diff.identical, diff.report()
    assert "arrays_differ" in diff.diff
    assert "root['roman']['data']" in diff.diff["arrays_differ"]


@pytest.mark.parametrize("modification", [None, "dtype", "names", "values", "meta"])
def test_compare_asdf_tables(tmp_path, modification):
    fn0 = tmp_path / "test0.asdf"
    fn1 = tmp_path / "test1.asdf"
    names = ("a",)
    values = [1, 2, 3]
    dtype = np.uint8
    t0 = astropy.table.Table([np.array(values, dtype=dtype)], names=names)
    if modification == "names":
        names = ("b",)
    if modification == "values":
        values = [4, 5, 6]
    if modification == "dtype":
        dtype = np.float64
    t1 = astropy.table.Table([np.array(values, dtype=dtype)], names=names)
    if modification == "meta":
        t1.meta["modified"] = True
    asdf.AsdfFile({"t": t0}).write_to(fn0)
    asdf.AsdfFile({"t": t1}).write_to(fn1)
    diff = compare_asdf(fn0, fn1)
    if modification is None:
        assert diff.identical, diff.report()
    else:
        assert not diff.identical, diff.report()
        assert "tables_differ" in diff.diff
