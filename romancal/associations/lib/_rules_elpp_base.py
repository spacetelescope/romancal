"""Base classes which define the ELPP Associations"""

from __future__ import annotations

import logging

from romancal.associations import libpath
from romancal.associations._association import _Association
from romancal.associations._exceptions import AssociationNotValidError
from romancal.associations._registry import RegistryMarker
from romancal.associations.lib._dms_base import DMSBaseMixin

__all__ = [
    "ASN_SCHEMA",
    "DMS_ELPP_Base",
]
# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
ASN_SCHEMA = RegistryMarker.schema(libpath("asn_schema_jw_level3.json"))


class DMS_ELPP_Base(DMSBaseMixin, _Association):
    """Basic class for DMS Level associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA.schema

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Other presumptions on the association
        if "constraints" not in self.data:
            self.data["constraints"] = "No constraints"
        if "asn_type" not in self.data:
            self.data["asn_type"] = "user_built"
        if "asn_id" not in self.data:
            self.data["asn_id"] = "a3001"
        if "target" not in self.data:
            self.data["target"] = "none"
        if "asn_pool" not in self.data:
            self.data["asn_pool"] = "none"
        if "skycell_wcs_info" not in self.data:
            self.data["skycell_wcs_info"] = "none"

    @property
    def current_product(self):
        return self.data["products"][-1]

    def __eq__(self, other):
        """Compare equality of two associations"""
        if isinstance(other, DMS_ELPP_Base):
            result = self.data["asn_type"] == other.data["asn_type"]
            result = result and (self.member_ids == other.member_ids)
            return result

        return NotImplemented

    def __ne__(self, other):
        """Compare inequality of two associations"""
        result = self.__eq__(other)
        if result is not NotImplemented:
            result = not result
        return result

    def _add_items(self, items, product_name=None, with_exptype=False, **kwargs):
        """Force adding items to the association

        Parameters
        ----------
        items : [object[, ...]]
            A list of items to make members of the association.

        product_name : str or None
            The name of the product to add the items to.
            If the product does not already exist, it will be created.
            If None, the default DMS ELPP naming
            conventions will be attempted.

        with_exptype : bool
            If True, each item is expected to be a 2-tuple with
            the first element being the item to add as `expname`
            and the second items is the `exptype`

        kwargs : dict
            Allows other keyword arguments used by other subclasses.

        Notes
        -----
        This is a low-level shortcut into adding members, such as file names,
        to an association. All defined shortcuts and other initializations are
        by-passed, resulting in a potentially unusable association.
        """
        if product_name is None:
            raise AssociationNotValidError("Product name needs to be specified")

        self.new_product(product_name)
        members = self.current_product["members"]
        for item in items:
            exptype = "science"
            if with_exptype:
                item, exptype = item
            member = {"expname": item, "exptype": exptype}
            members.append(member)
