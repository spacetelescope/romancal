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


from ._association import *
from ._association_io import *
from ._exceptions import *
from ._generate import *
from .lib.process_list import *
from ._load_asn import load_asn
from ._main import *
from ._pool import *
from ._registry import *
from .skycell_asn import run_skycell_asn, FileRecord, _cli, parse_visitID, _create_groups, _extract_visit_id, _fetch_filter_for, _save_association, _create_intersecting_skycell_index, _group_files_by_filter_for_skycell, _create_metadata, asn_from_list
from .asn_from_list import *
from .multiband_asn import *
