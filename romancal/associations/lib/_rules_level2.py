"""Association Definitions: DMS Level2b product associations"""

import logging

from romancal.associations._registry import RegistryMarker
from romancal.associations.lib._rules_elpp_base import *

__all__ = [
    "AsnMinxin_Lv2FOV",
    "AsnMixin_Lv2Image",
    "Asn_Lv2FOV",
    "Asn_Lv2GBTDSFull",
    "Asn_Lv2GBTDSPass",
    "Asn_Lv2Image",
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# flake8: noqa: F403, F405
# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Lv2FOV(AsnMixin_Lv2FOV, DMS_ELPP_Base):
    """Level2b Non-TSO Science Image Association

    Characteristics:
        - Association type: ``FOV``
        - Pipeline: ``mosaic``
        - Image-based science exposures
        - Science exposures for all 18 detectors
    """


@RegistryMarker.rule
class Asn_Lv2Image(AsnMixin_Lv2Image, DMS_ELPP_Base):
    """Level2b Non-TSO Science Image Association

    Characteristics:
        - Association type: ``image``
        - Pipeline: ``ELPP``
        - Image-based science exposures
        - Single science exposure
    """


@RegistryMarker.rule
class Asn_Lv2GBTDSPass(AsnMixin_Lv2GBTDSpass, DMS_ELPP_Base):
    """Level2 GBTDS Pass Science Image Association

    Characteristics:
        - Association type: ``image``
        - Pipeline: ``romancal_image``
        - Image-based science exposures
        - Collect all exposures in a pass
        - Non-TSO
    """


@RegistryMarker.rule
class Asn_Lv2GBTDSFull(AsnMixin_Lv2GBTDSfull, DMS_ELPP_Base):
    """Level2 GBTDS Full Science Image Association

    Characteristics:
        - Association type: ``image``
        - Pipeline: ``romancal_image``
        - Image-based science exposures
        - Collect all exposures in a season
        - Non-TSO
    """

    def get_exposure_type(self, item, default=None):
        """overrides super method to return `item["exp_type"]`"""

        return item["exp_type"]
