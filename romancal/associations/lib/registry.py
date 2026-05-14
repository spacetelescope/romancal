# Migrated from romancal/associations/registry.py
"""Association Registry (migrated)"""

import logging

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
