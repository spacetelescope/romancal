from contextlib import chdir

import pytest


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
