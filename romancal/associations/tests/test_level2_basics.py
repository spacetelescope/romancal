"""Test basic usage of Level2 associations"""

from romancal.associations import _generate
from romancal.associations.tests.helpers import (
    combine_pools,
    registry_level2_only,
    t_path,
)

REGEX_LEVEL2 = r"(?P<path>.+)(?P<type>_cal?)(?P<extension>\..+)"


def generate_from_pool(pool_path):
    """Generate associations from pools"""
    rules = registry_level2_only()
    pool = combine_pools(t_path(pool_path))
    asns = _generate(pool, rules)
    return asns


def test_level2_productname():
    asns = generate_from_pool("data/pool_002_wfi_image.csv")
    for asn in asns:
        for product in asn["products"]:
            science = [
                member
                for member in product["members"]
                if member["exptype"] == "science" or member["exptype"] == "wfi_image"
            ]
            if asn["asn_rule"] == "Asn_Lv2Image":
                assert len(science) == 2
            if asn["asn_rule"] == "Asn_Lv2FOV":
                assert len(science) == 18
