"""Test basic utils bytes2human"""

import pytest

from romancal.lib.basic_utils import bytes2human


test_data = [(1000, '1000B'),
             (1024, '1.0K'),
             (1024*10, '10.0K'),
             (100001221, '95.4M')]


@pytest.mark.parametrize("input_data, result", test_data)
def test_bytes2human(input_data, result):
    """ test the basic conversion """

    assert bytes2human(input_data) == result
