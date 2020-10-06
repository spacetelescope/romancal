# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
from pkg_resources import get_distribution, DistributionNotFound
if sys.version_info < (3, 6):
    raise ImportError("Romancal supports Python versions 3.6 and above.")  # pragma: no cover

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:  # pragma: no cover
    # package is not installed
    pass  # pragma: no cover
