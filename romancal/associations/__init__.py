"""Association Generation code

Code for generating associations, lists of association members and
relevant metadata describing a collection of exposures.

For more, see the :ref:`documentation overview <asn-overview>`.

"""

# flake8: noqa: F401

from ._association import _Association
from ._exceptions import AssociationNotValidError
from ._load_asn import load_asn
