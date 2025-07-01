#!/usr/bin/env python
"""
This test file is a bit different in that pytest should
not be imported within this file. It is designed to check
that all submodules are importable (without silently
depending on pytest or it's dependencies) and is run
in a separate tox environment where pytest is not installed.
"""

import importlib
import pkgutil

import romancal


def dependencies(package, exclude: [str]):
    return [
        module[1]
        for module in pkgutil.walk_packages(
            package.__path__, prefix=package.__name__ + "."
        )
        if not any(exclude_module in module[1] for exclude_module in exclude)
    ]


# dqflags is deprecated
MODULES = dependencies(romancal, exclude=["test", "time", "static_preview", "dqflags", "plot"])


def test_modules_import():
    for module in MODULES:
        importlib.import_module(module)


if __name__ == "__main__":
    test_modules_import()
