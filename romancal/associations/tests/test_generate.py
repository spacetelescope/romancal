"""Test basic generate operations"""

from romancal.associations import (
    AssociationPool,
    AssociationRegistry,
    generate,
    load_asn,
)
from romancal.associations.tests.helpers import t_path


def test_simple():
    """Test generate on simple registry"""
    registry = AssociationRegistry(
        [t_path("data/rules_basic.py")], include_default=False
    )
    pool = AssociationPool()
    pool["value"] = ["row1", "row2"]

    asns = generate(pool, registry)
    assert len(asns) == 1
    assert len(asns[0]["members"]) == 2


def test_unserialize():
    """Test basic unserializing"""
    asn_file = t_path("data/asn_mosaic.json")
    with open(asn_file) as asn_fp:
        asn = load_asn(asn_fp)
    assert isinstance(asn, dict)
