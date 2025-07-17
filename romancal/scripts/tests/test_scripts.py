import shutil

import pytest

SCRIPTS = [
    "roman_pointing_summary",
    "roman_set_telescope_pointing",
    "roman_v1_calculate",
    "verify_install_requires",
]


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_installed(script):
    assert shutil.which(script) is not None, f"`{script}` not installed"
