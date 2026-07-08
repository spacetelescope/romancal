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


def combine_pools(pools, **pool_kwargs):
    """Combine pools into a single pool

    Parameters
    ----------
    pools: str, astropy.table.Table, [str|Table, ...]
        The pools to combine. Either a singleton is
        passed or and iterable can be passed.
        The entries themselves can be either a file path
        or an astropy.table.Table-like object.

    pool_kwargs: dict
        Other keyword arguments to pass to AssociationPool.read

    Returns
    -------
    AssociationPool|astropy.table.Table
        The combined pool
    """
    if not is_iterable(pools):
        pools = [pools]
    just_pools = []
    for pool in pools:
        if not isinstance(pool, Table):
            pool = _AssociationPool.read(pool, **pool_kwargs)
        just_pools.append(pool)
    if len(just_pools) > 1:
        mega_pool = vstack(just_pools, metadata_conflicts="silent")
    else:
        mega_pool = just_pools[0].copy(copy_data=True)

    # Replace OBS_NUM and ASN_CANDIDATE_ID with actual numbers, if
    # necessary
    expnum = Counter(start=0)
    obsnum = Counter(start=0)
    acid = Counter(start=999)
    local_env = locals()
    global_env = globals()
    for row in mega_pool:
        mega_pool[row.index] = [
            parse_value(v, global_env=global_env, local_env=local_env) for v in row
        ]

    return mega_pool


def parse_value(v, global_env=None, local_env=None):
    """Evaluate if indicated"""
    if global_env is None:
        global_env = globals()
    if local_env is None:
        local_env = locals()

    result = v
    try:
        m = re.match("@!(.+)", v)
    except TypeError:
        pass
    else:
        if m:
            result = eval(m.group(1), global_env, local_env)
    return result


@contextmanager
def mkstemp_pool_file(pools, **pool_kwargs):
    """Make an actual pool file"""
    pool = combine_pools(pools, **pool_kwargs)
    with TemporaryDirectory() as path:
        pool_path = os.path.join(path, "pool")
        pool.write(
            pool_path,
            format="ascii",
            delimiter="|",
        )
        yield pool_path


def level2_rule_path():
    """Return the path to the level 2 rules"""
    return t_path("../lib/_rules_level2.py")


def registry_level2_only(global_constraints=None):
    """Get registry with only Level2 rules"""
    return _AssociationRegistry(
        definition_files=[level2_rule_path()],
        include_default=False,
        global_constraints=global_constraints,
    )
