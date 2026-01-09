import pytest


@pytest.fixture(scope="session")
def old_build_path():
    """
    Artifactory path of a previous build.

    The files for this build will be used to test
    pipeline backwards compatibility.
    """
    return "build/26Q1_B20"
