"""Test utility update_path"""

from romancal.associations.asn_from_list import asn_from_list
from romancal.associations.lib.rules_elpp_base import DMS_ELPP_Base
from romancal.associations.lib.update_path import update_path


def test_update_path_level2():
    members = ["a", "b", "c"]
    new_path = "new_path"
    product_name = "new_product"
    asn = asn_from_list(members, product_name=product_name, rule=DMS_ELPP_Base)
    update_path(asn, new_path)
    for product in asn["products"]:
        for member in product["members"]:
            assert member["expname"].startswith(new_path)


def test_update_path_level3():
    members = ["a", "b", "c"]
    new_path = "new_path"
    asn = asn_from_list(members, product_name="test")
    update_path(asn, new_path)
    for product in asn["products"]:
        for member in product["members"]:
            assert member["expname"].startswith(new_path)
