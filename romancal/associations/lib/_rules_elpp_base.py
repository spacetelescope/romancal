"""Base classes which define the ELPP Associations"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from romancal.associations import libpath
from romancal.associations._association import _Association
from romancal.associations._exceptions import AssociationNotValidError
from romancal.associations._registry import RegistryMarker
from romancal.associations.lib._constraint import Constraint, SimpleConstraint
from romancal.associations.lib._dms_base import (
    _EMPTY,
    IMAGE2_NONSCIENCE_EXP_TYPES,
    IMAGE2_SCIENCE_EXP_TYPES,
    SPEC2_SCIENCE_EXP_TYPES,
    DMSAttrConstraint,
    DMSBaseMixin,
)
from romancal.associations.lib._keyvalue_registry import KeyValueRegistryError
from romancal.associations.lib._member import Member
from romancal.associations.lib._process_list import ProcessList

if TYPE_CHECKING:
    pass

__all__ = [
    "ASN_SCHEMA",
    "AsnMixin_Lv2FOV",
    "AsnMixin_Lv2GBTDSfull",
    "AsnMixin_Lv2GBTDSpass",
    "AsnMixin_Lv2Image",
    "AsnMixin_Science",
    "AsnMixin_Spectrum",
    "Constraint",
    "Constraint_Base",
    "Constraint_Category",
    "Constraint_Expos",
    "Constraint_Filename",
    "Constraint_Image",
    "Constraint_Image_Science",
    "Constraint_Instrument",
    "Constraint_Obsnum",
    "Constraint_Optical_Path",
    "Constraint_Pass",
    "Constraint_Sequence",
    "Constraint_Single_Science",
    "Constraint_Spectral",
    "Constraint_Spectral_Science",
    "Constraint_SubCategory",
    "Constraint_Target",
    "Constraint_Tile",
    "DMSAttrConstraint",
    "DMS_ELPP_Base",
    "ProcessList",
    "SimpleConstraint",
]
# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The schema that these associations must adhere to.
ASN_SCHEMA = RegistryMarker.schema(libpath("asn_schema_jw_level3.json"))

# DMS file name templates
_LEVEL1B_REGEX = r"(?P<path>.+)(?P<type>_uncal)(?P<extension>\..+)"
_DMS_POOLNAME_REGEX = r"jw(\d{5})_(\d{8}[Tt]\d{6})_pool"

# Product name regex's
_REGEX_ACID_VALUE = r"(o\d{3}|(c|a)\d{4})"

# Exposures that should have received Level2b processing
LEVEL2B_EXPTYPES = []
LEVEL2B_EXPTYPES.extend(IMAGE2_SCIENCE_EXP_TYPES)
LEVEL2B_EXPTYPES.extend(IMAGE2_NONSCIENCE_EXP_TYPES)
LEVEL2B_EXPTYPES.extend(SPEC2_SCIENCE_EXP_TYPES)

# Association Candidates that should never make Level3 associations
INVALID_AC_TYPES = ["background"]


class DMS_ELPP_Base(DMSBaseMixin, _Association):
    """Basic class for DMS Level associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA.schema

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

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
            member = Member({"expname": item, "exptype": exptype}, item=item)
            members.append(member)

    def __repr__(self):
        try:
            _, json_repr = self.ioregistry["json"].dump(self)
        except KeyValueRegistryError:
            return str(self.__class__)
        return json_repr

    def __str__(self):
        result_list = []
        result_list.append(
            f"{self.asn_name} with {len(self.data['products'])} products"
        )
        result_list.append(f"Rule={self.data['asn_rule']}")
        result_list.append(self.data["constraints"])
        result_list.append("Products:")
        for product in self.data["products"]:
            result_list.append(
                f"\t{product['name']} with {len(product['members'])} members"
            )
        result = "\n".join(result_list)
        return result


# -----------------
# Basic constraints
# -----------------
class Constraint_Base(Constraint):
    """Select on program and instrument"""

    def __init__(self):
        super().__init__(
            [
                DMSAttrConstraint(
                    name="program",
                    sources=["program"],
                ),
                DMSAttrConstraint(
                    name="instrument",
                    sources=["instrume"],
                ),
            ],
            name="base",
        )


class Constraint_Instrument(Constraint):
    """Select on instrument"""

    def __init__(self):
        super().__init__(
            [
                DMSAttrConstraint(
                    name="instrument",
                    sources=["instrume"],
                ),
            ],
        )


class Constraint_Filename(DMSAttrConstraint):
    """Select on visit number"""

    def __init__(self):
        super().__init__(
            name="Filename",
            sources=["filename"],
        )


class Constraint_Expos(DMSAttrConstraint):
    """Select on exposure number"""

    def __init__(self):
        super().__init__(
            name="exposure_number",
            sources=["nexpsur"],
            force_unique=True,
            required=True,
        )


class Constraint_Tile(DMSAttrConstraint):
    """Select on exposure number"""

    def __init__(self):
        super().__init__(
            name="tile",
            sources=["tile"],
            force_unique=True,
            required=True,
        )


class Constraint_Image(DMSAttrConstraint):
    """Select on exposure type"""

    def __init__(self):
        super().__init__(
            name="exp_type",
            sources=["exp_type"],
            value="wfi_image|wfi_wfsc",
        )


class Constraint_Obsnum(DMSAttrConstraint):
    """Select on OBSNUM"""

    def __init__(self):
        super().__init__(
            name="obs_num",
            sources=["obs_num"],
            force_unique=False,
            required=False,
        )


class Constraint_Optical_Path(Constraint):
    """Select on optical path"""

    def __init__(self):
        super().__init__(
            [
                DMSAttrConstraint(
                    name="opt_elem",
                    sources=["opt_elem"],
                    required=False,
                )
            ]
        )


class Constraint_SubCategory(DMSAttrConstraint):
    """Select on SUBCAT (proposal subcategory) in the mock pool file"""

    def __init__(self):
        super().__init__(
            name="sub_cat",
            sources=["subcat"],
            force_unique=True,
            required=True,
        )


class Constraint_Category(DMSAttrConstraint):
    """Select on SUBCAT (proposal subcategory) in the mock pool file"""

    def __init__(self):
        super().__init__(
            name="category",
            sources=["cat"],
            force_unique=True,
            required=True,
        )


class Constraint_Pass(Constraint):
    """Select on pass number"""

    def __init__(self):
        super().__init__(
            [
                DMSAttrConstraint(
                    name="pass",
                    sources=["pass"],
                    required=True,
                    force_unique=True,
                )
            ],
            name="pass",
        )


class Constraint_Sequence(Constraint):
    """Select on pass number"""

    def __init__(self):
        super().__init__(
            [
                DMSAttrConstraint(
                    name="sequence",
                    sources=["sequence"],
                    required=True,
                    force_unique=True,
                )
            ],
            name="sequence",
        )


class Constraint_Image_Science(DMSAttrConstraint):
    """Select on science images"""

    def __init__(self):
        super().__init__(
            name="exp_type",
            sources=["exp_type"],
            value="|".join(IMAGE2_SCIENCE_EXP_TYPES),
        )


class Constraint_Single_Science(SimpleConstraint):
    """Allow only single science exposure

    Parameters
    ----------
    has_science_fn : func
        Function to determine whether the association
        has a science member already. No arguments are provided.

    sc_kwargs : dict
        Keyword arguments to pass to the parent class `SimpleConstraint`

    Notes
    -----
    The `has_science_fn` is further wrapped in a lambda function
    to provide a closure. Otherwise if the function is a bound method,
    that method may end up pointing to an instance that is not calling
    this constraint.
    """

    def __init__(self, has_science_fn, **sc_kwargs):
        super().__init__(
            name="single_science",
            value=False,
            sources=lambda item: has_science_fn(),
            **sc_kwargs,
        )


class Constraint_Spectral(DMSAttrConstraint):
    """Constrain on spectral exposure types"""

    def __init__(self):
        super().__init__(
            name="exp_type",
            sources=["exp_type"],
            value="wfi_grism|wfi_prism",
            force_unique=False,
        )


class Constraint_Spectral_Science(Constraint):
    """Select on spectral science

    Parameters
    exclude_exp_types : [exp_type[, ...]]
        List of exposure types to not consider from
        from the general list.
    """

    def __init__(self, exclude_exp_types=None):
        if exclude_exp_types is None:
            general_science = SPEC2_SCIENCE_EXP_TYPES
        else:
            general_science = set(SPEC2_SCIENCE_EXP_TYPES).symmetric_difference(
                exclude_exp_types
            )

        super().__init__(
            [
                DMSAttrConstraint(
                    name="exp_type",
                    sources=["exp_type"],
                    value="|".join(general_science),
                )
            ],
            reduce=Constraint.any,
        )


class Constraint_Target(DMSAttrConstraint):
    """Select on target

    Parameters
    ----------
    association: Association
        If specified, use the `get_exposure_type` method
        to as part of the target selection.
    """

    def __init__(self, association=None):
        if association is None:
            super().__init__(
                name="target",
                sources=["targname"],
            )
        else:
            super().__init__(
                name="target",
                sources=["targname"],
                onlyif=lambda item: association.get_exposure_type(item) != "background",
                force_reprocess=ProcessList.EXISTING,
                only_on_match=True,
            )


# -----------
# Base Mixins
# -----------


class AsnMixin_Science(DMS_ELPP_Base):
    """Basic science constraints"""

    def __init__(self, *args, **kwargs):
        # Setup target acquisition inclusion
        constraint_acqs = Constraint(
            [
                DMSAttrConstraint(
                    name="acq_obsnum",
                    sources=["obs_num"],
                    value=lambda: (
                        "(" + "|".join(self.constraints["obs_num"].found_values) + ")"
                    ),
                    force_unique=False,
                )
            ],
            name="acq_constraint",
            work_over=ProcessList.EXISTING,
        )

        # Put all constraints together.
        self.constraints = Constraint(
            [
                Constraint_Base(),
                DMSAttrConstraint(sources=["is_imprt"], force_undefined=True),
                Constraint(
                    [
                        Constraint(
                            [self.constraints, Constraint_Obsnum()], name="rule"
                        ),
                        constraint_acqs,
                    ],
                    name="acq_check",
                    reduce=Constraint.any,
                ),
            ],
            name="dmsbase_top",
        )

        super().__init__(*args, **kwargs)


class AsnMixin_Spectrum(AsnMixin_Science):
    """All things that are spectrum"""


# ---------------------------------------------
# Mixins to define the broad category of rules.
# ---------------------------------------------
class AsnMixin_Lv2FOV:
    """Level 2 Image association base"""


class AsnMixin_Lv2Image:
    """Level 2 Image association base"""


class AsnMixin_Lv2GBTDSpass:
    """Level 2 GBTDS association base"""


class AsnMixin_Lv2GBTDSfull:
    """Level 2 GBTDS association base"""
