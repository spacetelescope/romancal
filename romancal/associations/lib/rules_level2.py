"""Association Definitions: DMS Level2b product associations
"""
from collections import deque
import logging

from romancal.associations.exceptions import AssociationNotValidError
from romancal.associations.registry import RegistryMarker
from romancal.associations.lib.constraint import (Constraint, SimpleConstraint)
from romancal.associations.lib.dms_base import (
    format_list,
    item_getattr,
)
from romancal.associations.lib.member import Member
from romancal.associations.lib.process_list import ListCategory
from romancal.associations.lib.utilities import (getattr_from_list, getattr_from_list_nofail)
from romancal.associations.lib.rules_elpp_base import *
#from romancal.associations.lib.rules_level3_base import DMS_Level3_Base

__all__ = [
    'Asn_Lv2Image',
    'Asn_Lv2Spec',
    'Asn_Lv2WFSC',
    'AsnMixin_Lv2Image'
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------
# Start of the User-level rules
# --------------------------------
@RegistryMarker.rule
class Asn_Lv2Image(
        AsnMixin_Lv2Image,
        DMS_ELPP_Base
):
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
        self.constraints = Constraint([
            Constraint_Base(),
            Constraint_Target(),
            Constraint_Expos(),
            Constraint_Tile()
        ])

        # Now check and continue initialization.
        super(Asn_Lv2Image, self).__init__(*args, **kwargs)



    def get_exposure_type(self, item, default='science'):
        """Modify exposure type depending on dither pointing index

        Behaves as the superclass method. However, if the constraint
        `is_current_patt_num` is True, mark the exposure type as
        `background`.
        """
        #exp_type = super(Asn_Lv2MIRLRSFixedSlitNod, self).get_exposure_type(
        #    item, default
        #)
        #if exp_type == 'science' and self.constraints['is_current_patt_num'].matched:
        #    exp_type = 'background'

        exp_type = item['exp_type']

        if exp_type == 'wfi_image':
            exp_type == 'science'

        return exp_type

#    def _init_hook(self, item):
#        """Post-check and pre-add initialization"""
#
#        super(Asn_Lv2WFSC, self)._init_hook(item)
#        self.data['asn_type'] = 'wfs-image2'
