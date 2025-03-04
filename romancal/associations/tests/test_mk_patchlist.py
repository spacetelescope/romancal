"""Tests for mk_patchlist"""

import pytest

# from romancal.associations import Association, AssociationRegistry, load_asn
import romancal.associations.mk_patchlist as mk_patchlist
from romancal.associations.mk_patchlist import _cli


def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        _cli([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        _cli(["-o", "test_asn.json"])
