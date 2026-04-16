import pytest

from romancal.datamodels.migration import MigrationWarning


@pytest.fixture(scope="module", params=["build/26Q1_B20", "build/26Q2_B21"])
def old_rtdata_module(rtdata_module, request):
    rtdata_module._env = request.param
    with pytest.warns(MigrationWarning, match="hga_move"):
        yield rtdata_module
