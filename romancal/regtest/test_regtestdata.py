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
