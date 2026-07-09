"""Association Registry"""

import importlib.util
import logging
from inspect import getmembers, isclass, isfunction, ismethod, ismodule
from os.path import basename, expanduser, expandvars

__all__ = ["RegistryMarker"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class RegistryMarker:
    """Mark rules, and modules for inclusion into a registry"""

    class Schema:
        def __init__(self, obj):
            self._asnreg_role = "schema"
            self._asnreg_schema = obj
            RegistryMarker.mark(self)

        @property
        def schema(self):
            return self._asnreg_schema

    @staticmethod
    def mark(obj):
        """Mark that an object should be part of the registry

        Parameters
        ----------
        obj : object
            The object to mark

        Returns
        -------
        obj
            Object that has been marked. Returned to enable
            use as a decorator.

        Notes
        -----
        The following attributes are added to the object:

        - _asnreg_mark : True
              Attribute added to object and is set to True

        - _asnreg_role : str or None
              If not already assigned, the role is left
              unspecified using None.

        """
        obj._asnreg_marked = True
        obj._asnreg_role = getattr(obj, "_asnreg_role", None)
        return obj

    @staticmethod
    def rule(obj):
        """Mark object as rule

        Parameters
        ----------
        obj : object
            The object that should be treated as a rule

        Returns
        -------
        obj : object
            Return object to enable use as a decorator.

        Notes
        -----
        The following attributes are added to the object:

        - _asnreg_role : 'rule'
              Attributed added to object and set to `rule`

        - _asnreg_mark : True
              Attributed added to object and set to True
        """
        obj._asnreg_role = "rule"
        RegistryMarker.mark(obj)
        return obj

    @staticmethod
    def schema(filename):
        """Mark a file as a schema source"""
        schema = RegistryMarker.Schema(filename)
        return schema

    @staticmethod
    def is_marked(obj):
        """Has an objected been marked?"""
        return hasattr(obj, "_asnreg_marked")


# Utilities
def import_from_file(filename):
    """Import a file as a module

    Parameters
    ---------
    filename : str
        The file to import

    Returns
    -------
    module : python module
        The imported module
    """
    path = expandvars(expanduser(filename))
    module_name = basename(path).split(".")[0]
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def get_marked(module, predicate=None, include_bases=False):
    """Recursively get all executable objects

    Parameters
    ----------
    module : python module
        The module to examine

    predicate : bool func(object)
        Determinant of what gets returned.
        If None, all object types are examined

    include_bases : bool
        If True, include base classes not considered
        rules.

    Returns
    -------
    class object : generator
        A generator that will yield all class members in the module.
    """

    def is_method(obj):
        return isfunction(obj) or ismethod(obj)

    for name, obj in getmembers(module, predicate):
        if isclass(obj):
            yield from get_marked(obj, predicate=is_method)
            if RegistryMarker.is_marked(obj) or include_bases:
                yield name, obj
        elif RegistryMarker.is_marked(obj):
            if ismodule(obj):
                yield from get_marked(
                    obj, predicate=predicate, include_bases=include_bases
                )
            else:
                yield name, obj
