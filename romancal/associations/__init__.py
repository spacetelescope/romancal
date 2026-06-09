"""Association Generator

The Association Generator takes a list of items, an Association Pool, and
creates sub-lists of those items depending on each item's attributes. How the
sub-lists are created is defined by Association Rules.

For more, see the :ref:`documentation overview <asn-overview>`.

"""

# flake8: noqa: F401, F403

# Take version from the upstream package
from .. import __version__


# Utility
def libpath(filepath):
    """Return the full path to the module library."""
    from os.path import abspath, dirname, join

    return join(dirname(abspath(__file__)), "lib", filepath)


from ._association import _Association
from ._association_io import AssociationNotValidError, json
from ._generate import _generate
from ._load_asn import load_asn
from ._pool import _AssociationPool
from ._registry import _AssociationRegistry

