import pytest


@pytest.fixture(
    scope="module",
    params=["build/26Q1_B20", "build/26Q2_B21", "build/26Q3_B22.3", "build/26Q4_B23"],
)
def old_rtdata_module(rtdata_module, request):
    rtdata_module._env = request.param
    yield rtdata_module
