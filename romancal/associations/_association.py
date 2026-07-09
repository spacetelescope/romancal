"""Main module for associations"""

import json
import logging
from collections.abc import MutableMapping

import jsonschema

from . import __version__
from ._exceptions import AssociationNotValidError

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

    def __init__(
        self,
        version_id=None,
        target=None,
    ):
        self.version_id = version_id
        self.target = target
        self.data = {
            "asn_type": "None",
            "asn_rule": self.asn_rule,
            "version_id": self.version_id,
            "code_version": __version__,
            "target": self.target,
        }
        self.meta = {}

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

    def dump(self):
        """Serialize the association

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
            asn_filename = self.asn_name
            if not asn_filename.endswith(".json"):
                asn_filename = asn_filename + ".json"
            return (
                asn_filename,
                json.dumps(self.data, indent=4, separators=(",", ": ")),
            )
        raise AssociationNotValidError(f"Association {self} is not valid")

    @classmethod
    def load(cls, serialized, validate=True):
        """Marshall a previously serialized association

        Parameters
        ----------
        serialized : object
            The serialized form of the association.

        validate : bool
            Validate against the class' defined schema, if any.

        Returns
        -------
        association : dict
            The association.

        Raises
        ------
        AssociationNotValidError
            Cannot create or validate the association.

        Notes
        -----
        The `serialized` object can be either a string or
        a file object containing the string.
        """
        try:
            if isinstance(serialized, str):
                loader = json.loads
            else:
                # Presume a file object
                serialized.seek(0)
                loader = json.load
            asn = loader(serialized)
        except Exception as err:
            logger.debug(f'Error unserializing: "{err}"')
            raise AssociationNotValidError(
                f'Container is not JSON: "{serialized}"'
            ) from err

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
