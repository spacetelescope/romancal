"""Helpers for tests."""

# ruff: noqa
import os
import re
from collections import namedtuple
from contextlib import contextmanager
from glob import glob
from tempfile import TemporaryDirectory

import pytest
from astropy.table import Table, vstack

from romancal.associations import _AssociationPool, _AssociationRegistry, _generate
from romancal.associations.lib._counter import Counter
from romancal.associations.lib._utilities import is_iterable


def t_path(partial_path):
    """Construction the full path for test files"""
    test_dir = os.path.dirname(__file__)
    return os.path.join(test_dir, partial_path)
