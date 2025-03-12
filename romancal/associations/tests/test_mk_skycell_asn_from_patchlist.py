"""Tests for mk_patchlist"""

import pytest
import numpy as np
from types import NoneType

# from romancal.associations import Association, AssociationRegistry, load_asn
import romancal.associations.mk_skycell_asn_from_patchlist as mk_skycell_asn_from_patchlist
from romancal.associations.mk_skycell_asn_from_patchlist import _cli


def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        _cli([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        _cli(["-o", "test_asn.json"])


def override_patch_table(monkeypatch):
    """
    For the tests in this file, monkeypatch the global
    PATCH_TABLE to a smaller PATCH_SUBSET to allow these tests
    to run without access to the full patch table.
    """
    monkeypatch.setattr(pm, "PATCH_PATH_TABLE", PATCH_SUBSET)
    yield

@pytest.mark.parametrize(
    "input_value, output_type",
    [
        (np.int64(55), type(dict)),
        (55, NoneType),
        (55., NoneType),
        ('r000dp90x00y55', NoneType),
    ],
)

def test_get_projectioncell_wcs(input_value, output_type):
    """ Test for getting the projection information for wcs information"""

    output = type(mk_skycell_asn_from_patchlist.get_projectioncell_wcs(input_value))
    assert isinstance(output, type(output_type))
