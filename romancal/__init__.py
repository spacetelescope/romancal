# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pkg_resources import DistributionNotFound, get_distribution

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:  # pragma: no cover
    # package is not installed
    pass  # pragma: no cover
