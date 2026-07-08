"""
Define the I/O methods for Level 3 associations
"""

import json as json_lib
import logging

from ._association import _Association
from ._exceptions import AssociationNotValidError
from .lib._member import Member

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = []


# Define JSON encoder to convert `Member` to `dict`
class AssociationEncoder(json_lib.JSONEncoder):
    """Encode to handle Associations"""

    def default(self, obj):
        # Convert Member to a simple dict
        if isinstance(obj, Member):
            return obj.data


@_Association.ioregistry
class json:
    """Load and store associations as JSON"""

    @staticmethod
    def load(cls, serialized):
        """Unserialize an association from JSON

        Parameters
        ----------
        cls : class
            The class from which further information will be gathered
            and possibly instantiated.

        serialized : str or file object
            The JSON to read

        Returns
        -------
        association : dict
            The association

        Raises
        ------
        AssociationNotValidError
            Cannot create or validate the association.
        """
        if isinstance(serialized, str):
            loader = json_lib.loads
        else:
            # Presume a file object
            serialized.seek(0)
            loader = json_lib.load
        try:
            asn = loader(serialized)
        except Exception as err:
            logger.debug(f'Error unserializing: "{err}"')
            raise AssociationNotValidError(
                f'Container is not JSON: "{serialized}"'
            ) from err

        return asn

    @staticmethod
    def dump(asn):
        """Create JSON representation.

        Parameters
        ----------
        asn : Association
            The association to serialize

        Returns
        -------
        (name, str):
            Tuple where the first item is the suggested
            Name for the JSON file.
            Second item is the string containing the JSON serialization.
        """
        asn_filename = asn.asn_name
        if not asn_filename.endswith(".json"):
            asn_filename = asn_filename + ".json"
        return (
            asn_filename,
            json_lib.dumps(
                asn.data, cls=AssociationEncoder, indent=4, separators=(",", ": ")
            ),
        )
