import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "stpsf: mark test as requiring stpsf data to run"
    )


def pytest_addoption(parser):
    parser.addoption(
        "--stpsf", action="store_true", default=False, help="run stpsf tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--stpsf"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_stpsf = pytest.mark.skip(reason="need --stpsf option to run")
    for item in items:
        if "stpsf" in item.keywords:
            item.add_marker(skip_stpsf)
