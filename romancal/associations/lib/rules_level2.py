"""Association Definitions: DMS Level2b product associations"""

import logging

from romancal.associations.lib.constraint import Constraint
from romancal.associations.lib.rules_elpp_base import *
from romancal.associations.registry import RegistryMarker

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

    def __init__(self, *args, **kwargs):
        # Setup constraints
        self.constraints = Constraint(
            [
                Constraint_Base(),
                Constraint_Target(),
                Constraint_Filename(),
            ]
        )

        # Now check and continue initialization.
        super().__init__(*args, **kwargs)


@RegistryMarker.rule
class Asn_Lv2Image(AsnMixin_Lv2Image, DMS_ELPP_Base):
    """Level2b Non-TSO Science Image Association

    Characteristics:
        - Association type: ``image``
        - Pipeline: ``ELPP``
        - Image-based science exposures
        - Single science exposure
    """

    def __init__(self, *args, **kwargs):
        # Setup constraints
        self.constraints = Constraint(
            [
                Constraint_Base(),
                Constraint_Target(),
                Constraint(
                    [
                        Constraint_Expos(),
                    ],
                    reduce=Constraint.any,
                ),
                Constraint_Optical_Path(),
                Constraint_Sequence(),
                Constraint_Pass(),
                Constraint_Tile(),
            ]
        )

        # Now check and continue initialization.
        super().__init__(*args, **kwargs)


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

    def __init__(self, *args, **kwargs):
        # Setup constraints
        self.constraints = Constraint(
            [
                Constraint_Base(),
                Constraint_Target(),
                Constraint_Category(),
                Constraint_Pass(),
                Constraint_Optical_Path(),
                Constraint_SubCategory(),
            ]
        )

        # Now check and continue initialization.
        super().__init__(*args, **kwargs)


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

    def __init__(self, *args, **kwargs):
        # Setup constraints
        self.constraints = Constraint(
            [
                Constraint_Base(),
                Constraint_Target(),
                Constraint_Category(),
                Constraint_Optical_Path(),
                Constraint_SubCategory(),
            ]
        )

        # Now check and continue initialization.
        super().__init__(*args, **kwargs)

    def get_exposure_type(self, item, default=None):
        """overrides super method to return `item["exp_type"]`"""

        return item["exp_type"]
