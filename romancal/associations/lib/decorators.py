"""
Generic decorators

Notes
-----
Lifted from the `glue-viz`_ project.

.. _glue-viz: https://github.com/glue-viz/glue
"""

from functools import wraps

__all__ = ["singleton"]


def _make_key(args, kwargs):
    return args, frozenset(kwargs.items())

def singleton(cls):
    """Turn a class into a singleton, such that new objects
    in this class share the same instance"""
    instances = {}
