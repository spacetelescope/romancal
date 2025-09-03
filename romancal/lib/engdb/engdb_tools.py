"""
Access the Roman Engineering Mnemonic Database.

The engineering mnemonics are provided by multiple services,
all of which require a level of authentication.

For non-operational use, the providing service is through the MAST AUI website
(TBD).

Authorization can be requested through https://auth.mast.stsci.edu/ site.

**Interface**

The primary entry point is the function `romancal.lib.engdb_tools.engdb_service`.
This function returns a `romancal.lib.engdb_lib.EngdbABC` connection object. Using
this object, values for a mnemonic covering a specified time range can be
retrieved using the :meth:`~romancal.lib.engdb_lib.EngdbABC.get_values` method.

By default, only values inclusively between the time end points are returned.
Depending on the frequency a mnemonic is updated, there can be no values. If
values are always desired, the nearest, bracketing values outside the time
range can be requested.

List of available services can be retrieved using `engdb_tools.AVAILABLE_SERVICES`

The primary service is the "mast" service. See `engdb_mast` for details.

.. warning::

    Many mnemonics are updated very quickly, up to 16Hz. When in doubt, specify a
    very short time frame, and request bracketing values. Otherwise, the request
    can return a very large amount of data, risking timeout, unnecessary memory
    consumption, or access restrictions.

Examples
--------
The typical workflow is as follows:

.. code-block:: python

    from romancal.lib.engdb.engdb_tools import engdb_service

    service = engdb_service()  # By default, will use the public MAST service.

    values = service.get_values("sa_zattest2", "2021-05-22T00:00:00", "2021-05-22T00:00:01")
"""

import logging

from .engdb_edp import EngdbEDP
from .engdb_mast import EngdbMast

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Define the available services
AVAILABLE_SERVICES = {
    "edp": EngdbEDP,
    "mast": EngdbMast,
}

# Expected errors from service initialization
EXPECTED_ERRORS = (KeyError, RuntimeError, TypeError)

__all__ = ["engdb_service"]


def engdb_service(service=None, **service_kwargs):
    """
    Provide access to the Roman Engineering Database.

    Access can be either through the public MAST API
    or by direct connection to the database server.

    Parameters
    ----------
    service : str or None
        The specific service to use. If None, first working service will be used.

    **service_kwargs : dict
        Service-specific keyword arguments. Refer to the concrete implementations
        of `~romancal.lib.engdb_lib.EngdbABC`.

    Returns
    -------
    engdb : `~romancal.lib.engdb_lib.EngdbABC`
        The engineering database service to use.
    """
    if service:
        try:
            engdb = AVAILABLE_SERVICES[service](**service_kwargs)
        except EXPECTED_ERRORS as excp:
            raise RuntimeError(f"Service {service} instantiation failed") from excp
        return engdb

    # No service was specified. Try until one is found.
    for name, service_class in AVAILABLE_SERVICES.items():
        try:
            engdb = service_class(**service_kwargs)
        except EXPECTED_ERRORS as excp:
            logger.debug("Service %s is unavailable.", name)
            logger.debug("Exception: %s", excp)
        else:
            # Found the working service. Continue on.
            break
    else:
        raise RuntimeError("No active engineering database service can be found.")

    # Service is in hand.
    return engdb
