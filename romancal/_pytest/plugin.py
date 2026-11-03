"""
These are automatically available when ``ci_watson`` is used as a
pytest plugin.
"""

import os

import pytest

from romancal._pytest.resource_tracker import ResourceTracker

__all__ = []


def pytest_addoption(parser):
    """
    These pytest hooks allow us to mark tests and run the marked tests with
    specific command line options.
    """
    # Add option to use big data sets
    parser.addoption(
        "--bigdata", action="store_true", help="use big data sets (intranet)"
    )

    # Choose to test under dev or stable
    parser.addoption("--env", default="dev", help="specify what environment to test")

    # Data file input/output source/destination customization.
    parser.addini(
        "inputs_root",
        "Root dir (or data repository name) for test input files.",
        default=None,
    )

    parser.addini(
        "results_root",
        "Root dir (or data repository name) for test result/output files.",
        default=None,
    )


def pytest_configure(config):
    config.getini("markers").append("bigdata: Run tests that require intranet access")


def pytest_runtest_setup(item):
    if "bigdata" in item.keywords and not item.config.getvalue("bigdata"):
        pytest.skip("need --bigdata option to run")


@pytest.fixture(scope="function")
def _jail(tmp_path):
    """Perform test in a pristine temporary working directory."""
    old_dir = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield str(tmp_path)
    finally:
        os.chdir(old_dir)


@pytest.fixture(scope="module")
def resource_tracker():
    """Fixture to return the current module-scoped ResourceTracker.

    Use by calling ``track`` to generate a context in which resource
    usage will be tracked.

    .. code-block:: python

        with resource_tracker.track():
            # do stuff

    For resources used during tests providing a function-scoped
    request fixture result as the log argument will also log the
    used resources to the junit results.xml.

    .. code-block:: python

        def test_something(resource_tracker, request):
            with resource_tracker.track(log=request):
                # do stuff

    For resources used during fixtures the tracked resources
    can be logged in a separate test using ``log_tracked_resources``.
    """
    return ResourceTracker()


@pytest.fixture()
def log_tracked_resources(resource_tracker, request):
    """Fixture to log resources tracked by ``resource_tracker``.

    .. code-block:: python

        @pytest.fixture
        def my_fixture(resource_tracker):
            with resource_tracker.track():
                # do stuff

        def test_write_log(log_tracked_resources, my_fixture):
            log_tracked_resources()
    """

    def callback():
        resource_tracker.log(request)

    yield callback
