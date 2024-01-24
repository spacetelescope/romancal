import asdf
import astropy.table
import numpy as np
import pytest
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.assign_wcs.assign_wcs_step import load_wcs
from romancal.regtest.regtestdata import compare_asdf


def _add_wcs(tmp_path, model):
    dfn = tmp_path / "wcs_distortion.asdf"
    distortion_model = rdm.DistortionRefModel(maker_utils.mk_distortion())
    distortion_model.save(dfn)
    load_wcs(model, {"distortion": dfn})


@pytest.mark.parametrize("modification", [None, "small", "large"])
def test_compare_asdf(tmp_path, modification):
    # Need to use different directories for the files because the filenames are
    #    now always updated as part of writing a datamodel, see spacetelescope/datamodels#295
    file0 = tmp_path / "test0"
    file1 = tmp_path / "test1"
    file0.mkdir()
    file1.mkdir()

    file0 = file0 / "test.asdf"
    file1 = file1 / "test.asdf"

    l2 = rdm.ImageModel(maker_utils.mk_level2_image(shape=(100, 100)))
    _add_wcs(tmp_path, l2)
    l2.save(file0)
    atol = 0.0001
    if modification == "small":
        l2.data += (atol / 2) * l2.data.unit
    elif modification == "large":
        l2.data += (atol * 2) * l2.data.unit
    l2.save(file1)
    diff = compare_asdf(file0, file1, atol=atol)
    if modification == "large":
        assert not diff.identical, diff.report()
        assert "arrays_differ" in diff.diff
        assert "root['roman']['data']" in diff.diff["arrays_differ"]
    else:
        assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "modification",
    [
        None,
        "dtype",
        "names",
        "small_values",
        "large_values",
        "small_meta",
        "large_meta",
    ],
)
def test_compare_asdf_tables(tmp_path, modification):
    fn0 = tmp_path / "test0.asdf"
    fn1 = tmp_path / "test1.asdf"
    atol = 0.0001
    names = ("a",)
    dtype = np.float32
    values = np.array([1, 2, 3], dtype=dtype)
    t0 = astropy.table.Table([values], names=names)
    # use atol so value is a float
    t0.meta["value"] = atol
    if modification == "names":
        names = ("b",)
    if modification == "small_values":
        values = values + atol / 2
    if modification == "large_values":
        values = values + atol * 2
    if modification == "dtype":
        dtype = np.float64
    t1 = astropy.table.Table([np.array(values, dtype=dtype)], names=names)
    t1.meta["value"] = t0.meta["value"]
    if modification == "small_meta":
        t1.meta["value"] += atol / 2
    if modification == "large_meta":
        t1.meta["value"] += atol * 2
    asdf.AsdfFile({"t": t0}).write_to(fn0)
    asdf.AsdfFile({"t": t1}).write_to(fn1)
    diff = compare_asdf(fn0, fn1, atol=atol)
    if modification in (None, "small_values", "small_meta"):
        assert diff.identical, diff.report()
    else:
        assert not diff.identical, diff.report()
        assert "tables_differ" in diff.diff


def test_model_difference(tmp_path):
    fn0 = tmp_path / "a.asdf"
    fn1 = tmp_path / "b.asdf"
    ma = rdm.DistortionRefModel(maker_utils.mk_distortion())
    mb = rdm.LinearityRefModel(maker_utils.mk_linearity())
    ma.save(fn0)
    mb.save(fn1)
    diff = compare_asdf(fn0, fn1)
    assert not diff.identical
    assert (
        """'type_changes': {"root['roman']": {'new_type': <class 'roman_datamodels.stnode.DistortionRef'>"""
        in diff.report()
    )
