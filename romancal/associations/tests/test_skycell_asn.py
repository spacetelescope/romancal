"""Tests for skycell_asn"""

import pytest

import romancal.associations.skycell_asn as skycell_asn
from romancal.associations.skycell_asn import _cli


@pytest.mark.parametrize(
    "args",
    [
        [],
        ["-o", "test_asn.json"],
    ],
)
def test_cmdline_fails(args):
    """Exercise the command line interface"""
    with pytest.raises(SystemExit):
        _cli(args)


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


@pytest.mark.parametrize(
    "product_type,expected_key_count",
    [
        ("full", 1),
        ("visit", 2),
        ("pass", None),  # Don't check key count for "pass"
        (None, 1),
    ],
)
def test_create_groups_param(sample_filelist, product_type, expected_key_count):
    groups = skycell_asn._create_groups(sample_filelist, product_type)
    assert all(isinstance(v, list) for v in groups.values())
    if product_type in ("visit", "full", None):
        assert len(groups) == expected_key_count
    elif product_type == "pass":
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
        assert rec.filter_id.lower() == "f158"
        assert rec.skycell_indices == [1, 2]


@pytest.mark.parametrize(
    "fname,expected",
    [
        ("r0000101002003004005_0001_wfi10_cal.asdf", "0000101002003004005"),
        ("0000101002003004005_0001_wfi10_cal.asdf", "0000101002003004005"),
        ("0000101002003004005_0001_wfi10_cal.asdf", "0000101002003004005"),  # no 'r'
        ("r12345.asdf", "12345"),  # short filename
    ],
)
def test_extract_visit_id_variants(fname, expected):
    assert skycell_asn._extract_visit_id(fname) == expected


def test_fetch_filter_for():
    file_index = [
        skycell_asn.FileRecord("file1.asdf", [1, 2], "f158"),
        skycell_asn.FileRecord("file2.asdf", [2, 3], "f146"),
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

    meta = skycell_asn._create_metadata(
        ["file1.asdf"],
        "d1",
        "asnfile",
        "skycell1",
        {"foo": "bar"},
        "0000101002003004005",
    )
    assert meta["asn_type"] == "image"
    assert meta["program"] == "P1"
    assert meta["data_release_id"] == "d1"
    assert meta["target"] == "skycell1"
    assert meta["skycell_wcs_info"] == {"foo": "bar"}


@pytest.mark.parametrize(
    "skycell_index,expected",
    [
        (1, {"f158": {"file1.asdf", "file3.asdf"}, "f105": {"file5.asdf"}}),
        (3, {"f146": {"file2.asdf", "file4.asdf"}, "f105": {"file5.asdf"}}),
    ],
)
def test_group_files_by_filter_for_skycell(skycell_index, expected):
    file_list = [
        skycell_asn.FileRecord("file1.asdf", [1, 2], "f158"),
        skycell_asn.FileRecord("file2.asdf", [2, 3], "f146"),
        skycell_asn.FileRecord("file3.asdf", [1], "f158"),
        skycell_asn.FileRecord("file4.asdf", [3], "f146"),
        skycell_asn.FileRecord("file5.asdf", [1, 3], "f105"),
    ]
    result = skycell_asn._group_files_by_filter_for_skycell(file_list, skycell_index)
    for k in expected:
        assert set(result[k]) == expected[k]
    assert set(result.keys()) == set(expected.keys())


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
