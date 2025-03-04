"""Tests for mk_patchlist"""

import pytest
import numpy as np
from types import NoneType

# from romancal.associations import Association, AssociationRegistry, load_asn
import romancal.associations.mk_skycellasn as mk_skycellasn
from romancal.associations.mk_skycellasn import _cli


def test_cmdline_fails():
    """Exercise the command line interface"""

    # No arguments
    with pytest.raises(SystemExit):
        _cli([])

    # Only the association file argument
    with pytest.raises(SystemExit):
        _cli(["-o", "test_asn.json"])


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

    output = type(mk_skycellasn.get_projectioncell_wcs(input_value))
    assert isinstance(output, type(output_type))
