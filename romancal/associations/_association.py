"""Main module for associations"""

import json
import logging
from collections.abc import MutableMapping

import jsonschema

from . import __version__
from ._exceptions import AssociationNotValidError
from .lib._constraint import Constraint
from .lib._ioregistry import IORegistry

__all__ = ["_Association"]


# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class _Association(MutableMapping):
    """Association Base Class

    Parameters
    ----------
    version_id : str or None
        Version ID to use in the name of this association.
        If None, nothing is added.

    Raises
    ------
    AssociationError
        If an item doesn't match.

    Attributes
    ----------
    instance : dict-like
        The instance is the association data structure.
        See `data` below

    meta : dict
        Information about the association.

    data : dict
        The association. The format of this data structure
        is determined by the individual associations and, if
        defined, validated against their specified schema.

    schema_file : str
        The name of the output schema that an association
        must adhere to.
    """

    registry = None
    """Registry this rule has been placed in."""

    DEFAULT_FORCE_UNIQUE = False
    """Default whether to force constraints to use unique values."""

    DEFAULT_REQUIRE_CONSTRAINT = True
    """Default require that the constraint exists or otherwise
    can be explicitly checked.
    """

    DEFAULT_EVALUATE = False
    """Default do not evaluate input values"""

    GLOBAL_CONSTRAINT = None
    """Global constraints"""

    INVALID_VALUES = None
    """Attribute values that indicate the
    attribute is not specified.
    """

    ioregistry = IORegistry()
    """The association IO registry"""

    def __init__(
        self,
        version_id=None,
        target=None,
    ):
        self.data = {}
        self.meta = {}

        self.version_id = version_id
        self.target = target
        self.data.update(
            {
                "asn_type": "None",
                "asn_rule": self.asn_rule,
                "version_id": self.version_id,
                "code_version": __version__,
                "target": self.target,
            }
        )

        # Setup constraints
        # These may be predefined by a rule.
        try:
            constraints = self.constraints
        except AttributeError:
            constraints = Constraint()
        if self.GLOBAL_CONSTRAINT is not None:
            constraints.append(self.GLOBAL_CONSTRAINT.copy())
        self.constraints = constraints

    @classmethod
    def create(cls, item, version_id=None):
        """Create association if item belongs

        Parameters
        ----------
        item : dict
            The item to initialize the association with.

        version_id : str or None
            Version ID to use in the name of this association.
            If None, nothing is added.

        Returns
        -------
        (association, reprocess_list)
            2-tuple consisting of:
                - association or None: The association or, if the item does not
                  match this rule, None
                - [ProcessList[, ...]]: List of items to process again.
        """
        asn = cls(version_id=version_id)

        matches, reprocess = asn.add(item)
        if not matches:
            return None, reprocess
        return asn, reprocess

    @property
    def asn_name(self):
        """Suggest filename for the association"""
        return "unnamed_association"

    @classmethod
    def _asn_rule(cls):
        return cls.__name__

    @property
    def asn_rule(self):
        """Name of the rule"""
        return self._asn_rule()

    @classmethod
    def validate(cls, asn):
        """Validate an association against this rule

        Parameters
        ----------
        asn : Association or association-like
            The association structure to examine

        Returns
        -------
        valid : bool
            True if valid. Otherwise the `AssociationNotValidError` is raised

        Raises
        ------
        AssociationNotValidError
            If there is some reason validation failed.

        Notes
        -----
        The base method checks against the rule class' schema
        If the rule class does not define a schema, a warning is issued
        but the routine will return True.
        """
        if not hasattr(cls, "schema_file"):
            logger.warning(f"Cannot validate: {cls} has no schema. Presuming OK.")
            return True

        if isinstance(asn, cls):
            asn_data = asn.data
        else:
            asn_data = asn

        with open(cls.schema_file) as schema_file:
            asn_schema = json.load(schema_file)

        try:
            jsonschema.validate(asn_data, asn_schema)
        except (AttributeError, jsonschema.ValidationError) as err:
            logger.debug("Validation failed:")
            logger.debug("%s", err)
            raise AssociationNotValidError("Validation failed") from err
        return True

    def dump(self, format="json", **kwargs):
        """Serialize the association

        Parameters
        ----------
        format : str
            The format to use to dump the association into.

        kwargs : dict
            List of arguments to pass to the registered
            routines for the current association type.

        Returns
        -------
        (name, serialized):
            Tuple where the first item is the suggested
            base name for the file.
            Second item is the serialization.

        Raises
        ------
        AssociationError
            If the operation cannot be done

        AssociationNotValidError
            If the given association does not validate.
        """
        if self.is_valid:
            return self.ioregistry[format].dump(self, **kwargs)

        raise AssociationNotValidError(f"Association {self} is not valid")

    @classmethod
    def load(cls, serialized, format=None, validate=True, **kwargs):
        """Marshall a previously serialized association

        Parameters
        ----------
        serialized : object
            The serialized form of the association.

        format : str or None
            The format to force. If None, try all available.

        validate : bool
            Validate against the class' defined schema, if any.

        kwargs : dict
            Other arguments to pass to the `load` method

        Returns
        -------
        association : Association
            The association.

        Raises
        ------
        AssociationNotValidError
            Cannot create or validate the association.

        Notes
        -----
        The `serialized` object can be in any format
        supported by the registered I/O routines. For example, for
        `json` and `yaml` formats, the input can be either a string or
        a file object containing the string.
        """
        if format is None:
            formats = [
                format_func for format_name, format_func in cls.ioregistry.items()
            ]
        else:
            formats = [cls.ioregistry[format]]

        for format_func in formats:
            try:
                asn = format_func.load(cls, serialized, **kwargs)
            except AssociationNotValidError:
                continue
            else:
                break
        else:
            raise AssociationNotValidError(
                f'Cannot translate "{serialized}" to an association'
            )

        # Validate
        if validate:
            cls.validate(asn)

        return asn

    @property
    def is_valid(self):
        """Check if association is valid"""
        try:
            self.__class__.validate(self)
        except AssociationNotValidError:
            return False
        return True

    def _add(self, item):
        """Add a item, association-specific"""
        raise NotImplementedError(
            "Association._add must be implemented by a specific association rule."
        )

    def _add_items(self, items, **kwargs):
        """Force adding items to the association

        Parameters
        ----------
        items : [object[, ...]]
            A list of items to make members of the association.

        Notes
        -----
        This is a low-level shortcut into adding members, such as file names,
        to an association. All defined shortcuts and other initializations are
        by-passed, resulting in a potentially unusable association.
        """
        try:
            self["members"].update(items)
        except KeyError:
            self["members"] = items

    # #################################################
    # Methods required for implementing MutableMapping
    # #################################################
    def __getitem__(self, key):
        return self.data[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.data[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.data[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __keytransform__(self, key):
        return key

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()

    def values(self):
        return self.data.values()
