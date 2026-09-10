"""Test basic generate operations"""

import importlib.resources

import romancal.associations
from romancal.associations import load_asn


def test_unserialize():
    """Test basic unserializing"""
    with open(
        importlib.resources.files(romancal.associations)
        / "tests"
        / "data"
        / "asn_mosaic.json"
    ) as asn_fp:
        asn = load_asn(asn_fp)
    assert isinstance(asn, dict)
