"""Association Generator

The Association Generator takes a list of items, an Association Pool, and
creates sub-lists of those items depending on each item's attributes. How the
sub-lists are created is defined by Association Rules.

For more, see the :ref:`documentation overview <asn-overview>`.

"""

# flake8: noqa: F401

from ._association import _Association
from ._exceptions import AssociationNotValidError
from ._load_asn import load_asn
