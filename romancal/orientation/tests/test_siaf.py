# Test SIAF access
from pathlib import Path

import pytest

# pysiaf is not a required dependency. If not present, ignore all this.
pysiaf = pytest.importorskip("pysiaf")

import romancal.orientation._siaf as siaf_lib  # noqa: E402


@pytest.mark.parametrize(
    "basepath, n_apers",
    [
        (None, 22),
        pytest.param(
            Path(__file__).parent / "data" / "siaf",
            5,
            marks=pytest.mark.xfail(
                reason="released pysiaf does not support basepath feature"
            ),
        ),
    ],
)
def test_open_siaf(basepath, n_apers):
    """Test opening siaf"""
    siaf = siaf_lib.open_siaf(basepath=basepath)
    assert len(siaf.apertures) == n_apers
