"""Tests for asn_from_list"""

import pytest

from romancal.associations import Association, load_asn
from romancal.associations.asn_from_list import asn_from_list
from romancal.associations.exceptions import AssociationNotValidError


def test_base_association():
    """Create the simplest of associations"""
    items = ["a", "b", "c"]
    asn = asn_from_list(items, rule=Association)
    assert asn["asn_rule"] == "Association"
    assert asn["asn_type"] == "None"
    assert asn["members"] == items


def test_base_roundtrip():
    """Write/read created base association"""
    items = ["a", "b", "c"]
    asn = asn_from_list(items, rule=Association)
    _, serialized = asn.dump()
    reloaded = load_asn(serialized, registry=None)
    assert asn["asn_rule"] == reloaded["asn_rule"]
    assert asn["asn_type"] == reloaded["asn_type"]
    assert asn["members"] == reloaded["members"]


def test_default_simple():
    """Default Level3 association"""
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
        [(item, type_) for item, type_ in items.items()],
        product_name=product_name,
        with_exptype=True,
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
