"""Association Definitions: DMS Level2b product associations
"""
import logging

from romancal.associations.lib.constraint import Constraint
from romancal.associations.lib.rules_elpp_base import *
from romancal.associations.registry import RegistryMarker

__all__ = ["Asn_Lv2Image", "Asn_Lv2GBTDSPass", "Asn_Lv2GBTDSFull", "AsnMixin_Lv2Image"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# flake8: noqa: F403, F405
# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Lv2Image(AsnMixin_Lv2Image, DMS_ELPP_Base):
    """Level2b Non-TSO Science Image Association

    Characteristics:
        - Association type: ``image2``
        - Pipeline: ``calwebb_image2``
        - Image-based science exposures
        - Single science exposure
        - Non-TSO
    """

    def __init__(self, *args, **kwargs):

        # Setup constraints
        self.constraints = Constraint(
            [
                Constraint_Base(),
                Constraint_Target(),
                Constraint_Expos(),
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

    def get_exposure_type(self, item, default="science"):
        """Modify exposure type depending on dither pointing index

        Behaves as the superclass method. However, if the constraint
        `is_current_patt_num` is True, mark the exposure type as
        `background`.
        """
        exp_type = item["exp_type"]

        if exp_type == "wfi_image":
            exp_type == "science"

        return exp_type
