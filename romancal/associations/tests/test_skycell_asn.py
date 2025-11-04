"""Tests for skycell_asn"""

from typing import ClassVar

import pytest

import romancal.associations.skycell_asn as skycell_asn
from romancal.associations.skycell_asn import _cli


def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        _cli([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        _cli(["-o", "test_asn.json"])


def test_parse_visitID():
    filelist1 = [
        "r0000101002003004005_0001_wfi10_cal.asdf",
    ]
    visitid_parts = skycell_asn.parse_visitID(filelist1[0][1:20])
    assert visitid_parts["Program"] == "00001"
    assert visitid_parts["Execution"] == "01"
    assert visitid_parts["Pass"] == "002"  # noqa: S105
    assert visitid_parts["Segment"] == "003"
    assert visitid_parts["Observation"] == "004"
    assert visitid_parts["Visit"] == "005"


@pytest.fixture
def sample_filelist():
    return [
        "r0000101002003004005_0001_wfi10_cal.asdf",
        "r0000101002003004005_0002_wfi10_cal.asdf",
        "r0000101002003004006_0001_wfi10_cal.asdf",
    ]


def test_create_groups_full(sample_filelist):
    groups = skycell_asn._create_groups(sample_filelist, "full")
    assert "full" in groups
    assert set(groups["full"]) == set(sample_filelist)


def test_create_groups_visit(sample_filelist):
    groups = skycell_asn._create_groups(sample_filelist, "visit")
    # Should group by visit id (first 19 chars after 'r')
    assert all(isinstance(v, list) for v in groups.values())
    assert len(groups) == 2  # two unique visit IDs in sample_filelist


def test_create_groups_pass(sample_filelist):
    groups = skycell_asn._create_groups(sample_filelist, "pass")
    assert all(isinstance(v, list) for v in groups.values())
    # Should group by pass key (positions :10)
    assert set(groups.keys())  # keys should not be empty


def test_create_intersecting_skycell_index(monkeypatch, sample_filelist):
    class DummyCalFile:
        class Meta:
            class Instrument:
                optical_element = "F158"

            instrument = Instrument()

            class WCS:
                pass

            wcs = WCS()

        meta = Meta()

        def close(self):
            pass

    monkeypatch.setattr(skycell_asn.rdm, "open", lambda fname: DummyCalFile())
    monkeypatch.setattr(skycell_asn.sm, "find_skycell_matches", lambda wcs: [1, 2])
    file_index = skycell_asn._create_intersecting_skycell_index(sample_filelist)
    assert len(file_index) == len(sample_filelist)
    for rec in file_index:
        assert rec[2].lower() == "f158"
        assert rec[1] == [1, 2]


def test_extract_visit_id():
    fname = "r0000101002003004005_0001_wfi10_cal.asdf"
    assert skycell_asn._extract_visit_id(fname) == "0000101002003004005"
    fname2 = "0000101002003004005_0001_wfi10_cal.asdf"
    assert skycell_asn._extract_visit_id(fname2) == "0000101002003004005"


def test_fetch_filter_for():
    file_index = [
        ["file1.asdf", [1, 2], "f158"],
        ["file2.asdf", [2, 3], "f146"],
    ]
    assert skycell_asn._fetch_filter_for("file1.asdf", file_index) == "f158"
    assert skycell_asn._fetch_filter_for("file2.asdf", file_index) == "f146"
    assert skycell_asn._fetch_filter_for("notfound.asdf", file_index) == "unknown"


def test_save_association(tmp_path):
    fname = tmp_path / "test_asn"
    content = '{"asn_type": "image"}'
    skycell_asn._save_association(str(fname), content)
    out_file = tmp_path / "test_asn_asn.json"
    assert out_file.exists()
    assert out_file.read_text() == content


def test_create_metadata(monkeypatch):
    class DummyASN(dict):
        def dump(self, *args, **kwargs):
            return None, '{"asn_type": "image"}'

    dummy_asn_obj = DummyASN()
    monkeypatch.setattr(
        skycell_asn.asn_from_list,
        "asn_from_list",
        lambda members, **kwargs: dummy_asn_obj,
    )
    monkeypatch.setattr(skycell_asn, "parse_visitID", lambda vid: {"Program": "P1"})

    class DummySkyCell:
        name = "skycell1"
        # need to use ClassVar here to avoid instance variable warning
        # (i.e., mutable default - dict - assigned as class attribute)
        wcs_info: ClassVar[dict] = {"foo": "bar"}

    meta = skycell_asn._create_metadata(
        ["file1.asdf"], "d1", "asnfile", DummySkyCell(), "0000101002003004005"
    )
    assert meta["asn_type"] == "image"
    assert meta["program"] == "P1"
    assert meta["data_release_id"] == "d1"
    assert meta["target"] == "skycell1"
    assert meta["skycell_wcs_info"] == {"foo": "bar"}


def test_create_groups_default_type(sample_filelist):
    """Test _create_groups with default (None) product_type returns 'full' group."""
    groups = skycell_asn._create_groups(sample_filelist, None)
    assert "full" in groups
    assert set(groups["full"]) == set(sample_filelist)


def test_group_files_by_visit_and_pass():
    """Test _group_files_by_visit and _group_files_by_pass group files correctly."""
    files = [
        "r0000101002003004005_0001_wfi10_cal.asdf",
        "r0000101002003004005_0002_wfi10_cal.asdf",
        "r0000101003003004006_0001_wfi10_cal.asdf",
    ]
    visit_groups = skycell_asn._group_files_by_visit(files)
    assert len(visit_groups) == 2
    for group in visit_groups.values():
        assert isinstance(group, list)
    pass_groups = skycell_asn._group_files_by_pass(files)
    assert set(pass_groups.keys()) == {"0000101002", "0000101003"}


def test_extract_visit_id_no_r():
    """Test _extract_visit_id when filename does not start with 'r'."""
    fname = "0000101002003004005_0001_wfi10_cal.asdf"
    assert skycell_asn._extract_visit_id(fname) == "0000101002003004005"


def test_extract_visit_id_short():
    """Test _extract_visit_id with a short filename."""
    fname = "r12345.asdf"
    assert skycell_asn._extract_visit_id(fname) == "12345"


def test_cli_parsing(monkeypatch):
    """Test _cli parses arguments and calls skycell_asn with correct values."""
    called = {}

    def fake_skycell_asn(filelist, output_file_root, product_type, data_release_id):
        called.update(
            {
                "filelist": filelist,
                "output_file_root": output_file_root,
                "product_type": product_type,
                "data_release_id": data_release_id,
            }
        )

    monkeypatch.setattr(skycell_asn, "skycell_asn", fake_skycell_asn)
    _cli_args = [
        "-o",
        "root",
        "--product-type",
        "visit",
        "--data-release-id",
        "d1",
        "f1.asdf",
        "f2.asdf",
    ]
    skycell_asn._cli(_cli_args)
    assert called["filelist"] == ["f1.asdf", "f2.asdf"]
    assert called["output_file_root"] == "root"
    assert called["product_type"] == "visit"
    assert called["data_release_id"] == "d1"


def test_group_files_by_filter_for_skycell():
    """Test _group_files_by_filter_for_skycell groups files by filter for a given skycell index."""
    # file_list: [filename, [skycell_indices], filter_id]
    file_list = [
        ["file1.asdf", [1, 2], "f158"],
        ["file2.asdf", [2, 3], "f146"],
        ["file3.asdf", [1], "f158"],
        ["file4.asdf", [3], "f146"],
        ["file5.asdf", [1, 3], "f105"],
    ]
    # Test for skycell_index = 1
    result = skycell_asn._group_files_by_filter_for_skycell(file_list, 1)
    assert set(result.keys()) == {"f158", "f105"}
    assert set(result["f158"]) == {"file1.asdf", "file3.asdf"}
    assert set(result["f105"]) == {"file5.asdf"}

    # Test for skycell_index = 3
    result = skycell_asn._group_files_by_filter_for_skycell(file_list, 3)
    assert set(result.keys()) == {"f146", "f105"}
    assert set(result["f146"]) == {"file2.asdf", "file4.asdf"}
    assert set(result["f105"]) == {"file5.asdf"}
