"""Test basic generate operations"""

from romancal.associations import load_asn
from romancal.associations.tests.helpers import t_path


def test_unserialize():
    """Test basic unserializing"""
    asn_file = t_path("data/asn_mosaic.json")
    with open(asn_file) as asn_fp:
        asn = load_asn(asn_fp)
    assert isinstance(asn, dict)
