"""Tests for asn_from_list"""

import pytest

from romancal.associations import _Association, load_asn
from romancal.associations._exceptions import AssociationNotValidError
from romancal.associations.asn_from_list import _cli, asn_from_list


def test_base_association():
    """Create the simplest of associations"""
    items = ["a", "b", "c"]
    with pytest.warns(DeprecationWarning, match="rule argument is deprecated"):
        asn = asn_from_list(items, rule=_Association)
    assert asn["asn_rule"] == "_Association"
    assert asn["asn_type"] == "None"
    assert asn["members"] == items


def test_base_roundtrip():
    """Write/read created base association"""
    items = ["a", "b", "c"]
    with pytest.warns(DeprecationWarning, match="rule argument is deprecated"):
        asn = asn_from_list(items, rule=_Association)
    _, serialized = asn.dump()
    reloaded = load_asn(serialized)
    assert asn["asn_rule"] == reloaded["asn_rule"]
    assert asn["asn_type"] == reloaded["asn_type"]
    assert asn["members"] == reloaded["members"]


def test_association_target():
    """Create the simple associations with target set"""
    items = ["a", "b", "c"]
    target_name = "270p65x48y69"
    product_name = "l3_target"
    rule_name = "_Association"
    with pytest.warns(DeprecationWarning, match="rule argument is deprecated"):
        asn = asn_from_list(
            items,
            rule=_Association,
            product_name=product_name,
            target=target_name,
            version_id="c55",
        )
    assert asn["asn_rule"] == rule_name
    assert asn["asn_type"] == "None"
    assert asn["members"] == items
    assert asn["target"] == target_name


def test_default_simple():
    """Default ELPP association"""
    product_name = "test_product"
    items = ["a", "b", "c"]
    asn = asn_from_list(items, product_name=product_name)
    assert asn["asn_rule"] == "DMS_ELPP_Base"
    assert asn["asn_type"] == "None"
    assert len(asn["products"]) == 1
    product = asn["products"][0]
    assert product["name"] == product_name
    assert len(product["members"]) == len(items)
    for member in product["members"]:
        assert member["expname"] in items
        assert member["exptype"] == "science"


def test_default_with_type():
    """ELPP association with types specified"""
    product_name = "test_product"
    items = {"a": "science", "b": "target_acq", "c": "somethingelse"}
    asn = asn_from_list(
        list(items.items()), product_name=product_name, with_exptype=True
    )
    assert asn["asn_rule"] == "DMS_ELPP_Base"
    assert asn["asn_type"] == "None"
    assert len(asn["products"]) == 1
    product = asn["products"][0]
    assert product["name"] == product_name
    assert len(product["members"]) == len(items)
    for member in product["members"]:
        assert member["expname"] in items
        assert member["exptype"] == items[member["expname"]]


def test_default_fail():
    """Test default DMS_ELPP_Base fail

    A product name needs to be included, but is not.
    """
    items = ["a"]
    with pytest.raises(AssociationNotValidError):
        _ = asn_from_list(items)


def test_default_roundtrip():
    """Create/Write/Read a ELPP association"""
    product_name = "test_product"
    items = {"a": "science", "b": "target_acq", "c": "somethingelse"}
    asn = asn_from_list(
        [(item, type_) for item, type_ in items.items()],
        product_name=product_name,
        with_exptype=True,
    )
    _, serialized = asn.dump()
    reloaded = load_asn(serialized)
    assert asn["asn_rule"] == reloaded["asn_rule"]
    assert asn["asn_type"] == reloaded["asn_type"]
    assert len(asn["products"]) == len(reloaded["products"])


def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        _cli([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        _cli(["-o", "test_asn.json"])


def test_cmdline_success(tmp_path):
    """Create ELPP associations in different formats"""
    path = tmp_path / "test_asn.json"
    product_name = "test_product"
    inlist = ["a", "b", "c"]
    args = ["-o", str(path), "--product-name", product_name]
    args = args + inlist
    return_code = _cli(args)
    with path.open() as fp:
        asn = load_asn(fp)
    assert len(asn["products"]) == 1
    assert asn["products"][0]["name"] == product_name
    members = asn["products"][0]["members"]
    expnames = [member["expname"] for member in members]
    assert inlist == expnames
    assert not return_code


def test_cmdline_change_rules(tmp_path):
    """Command line change the rule"""
    rule = "Asn_Lv2Image"
    path = tmp_path / "test_asn.json"
    inlist = ["a", "b", "c"]
    args = [
        "-o",
        str(path),
        "-r",
        rule,
        "--product-name",
        "test",
    ]
    args = args + inlist
    with pytest.warns(UserWarning, match="never supported"):
        _cli(args)
    with path.open() as fp:
        asn = load_asn(fp)
    # assert inlist == asn['members']
    assert inlist[0] == asn["products"][0]["members"][0]["expname"]
    assert inlist[1] == asn["products"][0]["members"][1]["expname"]
    assert inlist[2] == asn["products"][0]["members"][2]["expname"]


def test_api_list():
    """Test api call with simple list"""
    product_name = "test_product"
    inlist = ["a", "b", "c"]

    asn = asn_from_list(inlist, product_name=product_name)
    assert len(asn["products"]) == 1
    assert asn["products"][0]["name"] == product_name
    members = asn["products"][0]["members"]
    expnames = [member["expname"] for member in members]
    assert inlist == expnames


def test_api_with_type():
    """Test api call with type tuple"""
    product_name = "test_product"
    inlist = [("a", "science"), ("b", "psf"), ("c", "science")]

    asn = asn_from_list(inlist, product_name=product_name, with_exptype=True)
    assert len(asn["products"]) == 1
    assert asn["products"][0]["name"] == product_name
    members = asn["products"][0]["members"]
    members_dict = {member["expname"]: member["exptype"] for member in members}
    for name, type_ in inlist:
        assert members_dict[name] == type_


@pytest.mark.parametrize(
    "expname, catalog",
    [
        ("a_cal.asdf", "a_cat.parquet"),
        (
            "r0000101001001001001_0001_wfi01_f158_cal.asdf",
            "r0000101001001001001_0001_wfi01_f158_cat.parquet",
        ),
        # a directory on the member is kept on the catalog
        ("sub/dir/a_cal.asdf", "sub/dir/a_cat.parquet"),
        # the separator of the replaced suffix is preserved
        ("a-cal.asdf", "a-cat.parquet"),
        # an unknown suffix is appended to, not replaced
        ("a_unknown.asdf", "a_unknown_cat.parquet"),
    ],
)
def test_api_tweakreg_catalogs(expname, catalog):
    """Catalog names are derived from the member expnames"""
    asn = asn_from_list([expname], product_name="test", _tweakreg_catalogs=True)
    assert asn["products"][0]["members"][0]["tweakreg_catalog"] == catalog


def test_api_tweakreg_catalogs_off_by_default():
    """Members carry no catalog unless one was requested"""
    asn = asn_from_list(["a_cal.asdf"], product_name="test")
    assert "tweakreg_catalog" not in asn["products"][0]["members"][0]


def test_cmdline_tweakreg_catalogs(tmp_path):
    """Catalog names survive a trip through the association file"""
    path = tmp_path / "test_asn.json"
    inlist = ["a_cal.asdf", "b_cal.asdf"]
    args = [
        "-o",
        str(path),
        "--product-name",
        "test",
        "--_tweakreg-catalogs",
    ]
    _cli(args + inlist)

    with path.open() as fp:
        asn = load_asn(fp)
    members = asn["products"][0]["members"]
    catalogs = [member["tweakreg_catalog"] for member in members]
    assert catalogs == ["a_cat.parquet", "b_cat.parquet"]
