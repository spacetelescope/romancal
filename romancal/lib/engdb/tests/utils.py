# Testing utilities
import pytest


def assert_xfail(condition, reason='Unexpected database contents. Check state of database.'):
    """Instead of just failing, mark as expected fail"""
    if not condition:
        pytest.xfail(reason=reason)
