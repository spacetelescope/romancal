# Migrated from romancal/associations/registry.py
"""Association Registry (migrated)"""

import importlib.util
import logging
from inspect import getmembers, isclass, isfunction, ismethod, ismodule
from os.path import basename, expanduser, expandvars

from romancal.associations import libpath
from romancal.associations._exceptions import AssociationError, AssociationNotValidError
# If callback_registry is needed, migrate it as well
# from romancal.associations.lib.callback_registry import CallbackRegistry

__all__ = ["AssociationRegistry", "RegistryMarker"]

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Library files
_ASN_RULE = "association_rules.py"

class AssociationRegistry(dict):
    """The available associations (migrated stub)"""
    pass

class RegistryMarker:
    @staticmethod
    def schema(path):
        return path
