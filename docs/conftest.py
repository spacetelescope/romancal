import os
from contextlib import contextmanager

import pytest

# chdir contextmanager was added with Python 3.11
try:
    from contextlib import chdir
except ImportError:

    @contextmanager
    def chdir(path):
        old_path = os.getcwd()

        try:
            os.chdir(path)
            yield
        finally:
            os.chdir(old_path)


@pytest.fixture(autouse=True)
def _docdir(request):
    # Trigger ONLY for doctestplus.
    doctest_plugin = request.config.pluginmanager.getplugin("doctestplus")
    if isinstance(request.node.parent, doctest_plugin._doctest_textfile_item_cls):
        tmp_path = request.getfixturevalue("tmp_path")
        with chdir(tmp_path):
            yield
    else:
        yield
