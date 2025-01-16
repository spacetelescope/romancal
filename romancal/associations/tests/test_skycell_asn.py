"""Tests for skycell_asn"""

import pytest

# from romancal.associations import Association, AssociationRegistry, load_asn
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
