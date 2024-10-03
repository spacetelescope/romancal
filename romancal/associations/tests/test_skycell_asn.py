"""Tests for skycell_asn"""

import pytest

#from romancal.associations import Association, AssociationRegistry, load_asn
from romancal.associations.skycell_asn import Main

def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        Main([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        Main(["-o", "test_asn.json"])


