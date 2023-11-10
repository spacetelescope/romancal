"""Test basic utils"""

import numpy as np
import pytest
from astropy.table import Table

from romancal.lib.basic_utils import bytes2human, recarray_to_ndarray

test_data = [
    (1000, "1000B"),
    (1024, "1.0K"),
    (1024 * 10, "10.0K"),
    (100001221, "95.4M"),
]


@pytest.mark.parametrize("input_data, result", test_data)
def test_bytes2human(input_data, result):
    """test the basic conversion"""

    assert bytes2human(input_data) == result


def test_structured_array_utils():
    arrays = [np.arange(0, 10), np.arange(10, 20), np.arange(30, 40)]
    names = "a, b, c"

    recarr0 = np.core.records.fromarrays(
        arrays, names=names, formats=[arr.dtype for arr in arrays]
    )
    round_tripped = recarray_to_ndarray(recarr0)
    assert np.all(round_tripped == np.column_stack(arrays).astype("<f8"))

    astropy_table = Table(recarr0)
    assert np.all(
        np.array(arrays) == np.array([col.data for col in astropy_table.itercols()])
    )

    assert np.all(astropy_table.as_array() == recarr0)
