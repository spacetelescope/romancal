import shutil

import pytest

SCRIPTS = [
    "okify_regtests",
    "verify_install_requires",
]


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_installed(script):
    assert shutil.which(script) is not None, f"`{script}` not installed"
