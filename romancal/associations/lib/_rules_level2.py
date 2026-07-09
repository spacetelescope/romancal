"""Association Definitions: DMS Level2b product associations"""

import logging

from romancal.associations._registry import RegistryMarker
from romancal.associations.lib._rules_elpp_base import DMS_ELPP_Base

__all__ = [
    "Asn_Lv2Image",
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Lv2Image(DMS_ELPP_Base):
    """Level2b Non-TSO Science Image Association

    Characteristics:
        - Association type: ``image``
        - Pipeline: ``ELPP``
        - Image-based science exposures
        - Single science exposure
    """
