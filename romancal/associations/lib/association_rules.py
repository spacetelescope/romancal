"""Association Definitions: DMS-specific

Notes
-----
These associations are specifically defined for use in DMS.
"""
from romancal.associations.registry import RegistryMarker

from romancal.associations.lib import (rules_elpp_base)

RegistryMarker.mark(rules_elpp_base)
