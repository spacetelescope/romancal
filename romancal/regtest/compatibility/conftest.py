import pytest


@pytest.fixture(scope="module")
def old_rtdata_module(rtdata_module):
    rtdata_module._env = "build/26Q1_B20"
    yield rtdata_module
