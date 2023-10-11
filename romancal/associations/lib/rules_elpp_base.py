"""Base classes which define the ELPP Associations"""

import copy
import logging
import re
from collections import defaultdict
from os.path import basename, split, splitext

from stpipe.format_template import FormatTemplate

from romancal.associations import libpath
from romancal.associations.association import Association
from romancal.associations.exceptions import AssociationNotValidError
from romancal.associations.lib.acid import ACID
from romancal.associations.lib.constraint import Constraint, SimpleConstraint
from romancal.associations.lib.counter import Counter
from romancal.associations.lib.dms_base import (
    _EMPTY,
    IMAGE2_NONSCIENCE_EXP_TYPES,
    IMAGE2_SCIENCE_EXP_TYPES,
    SPEC2_SCIENCE_EXP_TYPES,
    WFI_DETECTORS,
    DMSAttrConstraint,
    DMSBaseMixin,
)
from romancal.associations.lib.keyvalue_registry import KeyValueRegistryError
from romancal.associations.lib.member import Member
from romancal.associations.lib.process_list import ProcessList
from romancal.associations.lib.product_utils import (
    prune_duplicate_associations,
    prune_duplicate_products,
)
from romancal.associations.lib.utilities import evaluate, is_iterable
from romancal.associations.registry import RegistryMarker

__all__ = [
    "ASN_SCHEMA",
    "AsnMixin_AuxData",
    "AsnMixin_Science",
    "AsnMixin_Spectrum",
    "AsnMixin_Lv2FOV",
    "AsnMixin_Lv2Image",
    "AsnMixin_Lv2GBTDSfull",
    "AsnMixin_Lv2GBTDSpass",
    "Constraint",
    "Constraint_Base",
    "Constraint_Category",
    "Constraint_Expos",
    "Constraint_Image",
    "Constraint_Instrument",
    "Constraint_Obsnum",
    "Constraint_Optical_Path",
    "Constraint_Pass",
    "Constraint_Spectral",
    "Constraint_SubCategory",
    "Constraint_Tile",
    "Constraint_Image_Science",
    "Constraint_Sequence",
    "Constraint_Single_Science",
    "Constraint_Spectral_Science",
    "Constraint_Target",
    "Constraint_Filename",
    "DMS_ELPP_Base",
    "DMSAttrConstraint",
    "ProcessList",
    "SimpleConstraint",
    "Utility",
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


class DMS_ELPP_Base(DMSBaseMixin, Association):
    """Basic class for DMS Level associations."""

    # Set the validation schema
    schema_file = ASN_SCHEMA.schema

    # Attribute values that are indicate the
    # attribute is not specified.
    INVALID_VALUES = _EMPTY

    # Make sequences type-dependent
    _sequences = defaultdict(Counter)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Initialize validity checks
        self.validity.update(
            {
                "has_science": {
                    "validated": True,
                    "check": lambda member: member["exptype"] == "science",
                },
            }
        )

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

    @property
    def current_product(self):
        return self.data["products"][-1]

    def members_by_type(self, member_type):
        """Get list of members by their exposure type"""
        member_type = member_type.lower()
        try:
            members = self.current_product["members"]
        except KeyError:
            result = []
        else:
            result = [
                member for member in members if member_type == member["exptype"].lower()
            ]

        return result

    def has_science(self):
        """Does association have a science member

        -------
        bool
            True if it does.
        """
        limit_reached = len(self.members_by_type("science")) >= 1
        return limit_reached

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

    @property
    def dms_product_name(self):
        """Define product name.

        Returns
        -------
        product_name : str
            The product name
        """
        return self._dms_product_name(self)

    @staticmethod
    def _dms_product_name(association):
        """Define product name.

        Parameters
        ----------
        association : `Association`
            Association to get the name from.

        Returns
        -------
        product_name : str
            The product name
        """
        target = association._get_target()

        instrument = association._get_instrument()

        opt_elem = association._get_opt_element()

        exposure = association._get_exposure()
        if exposure:
            exposure = "-" + exposure

        subarray = association._get_subarray()
        if subarray:
            subarray = "-" + subarray

        product_name = (
            "r{program}-{acid}" "_{target}" "_{instrument}" "_{opt_elem}{subarray}"
        )
        if "Full" in association.data["asn_rule"]:
            subarray = "Full"

        if "Pass" in association.data["asn_rule"]:
            subarray = "Pass"

        product_name = product_name.format(
            program=association.data["visit_id"],
            acid=association.acid.id,
            target=target,
            instrument=instrument,
            opt_elem=opt_elem,
            subarray=subarray,
            exposure=exposure,
        )

        return product_name.lower()

    def update_asn(self, item=None, member=None):
        """Update association meta information

        Parameters
        ----------
        item : dict or None
            Item to use as a source. If not given, item-specific
            information will be left unchanged.

        member : Member or None
            An association member to use as source.
            If not given, member-specific information will be update
            from current association/product membership.

        Notes
        -----
        If both `item` and `member` are given,
        information in `member` will take precedence.
        """
        super().update_asn(item=item, member=member)

        # Constraints
        self.data["constraints"] = str(self.constraints)

        # ID
        self.data["asn_id"] = self.acid.id

        # Target
        self.data["target"] = self._get_target()

        # Item-based information
        if item is not None:
            # Program
            if self.data["program"] == "noprogram":
                self.data["program"] = f"{item['program']:0>5s}"

            # Pool
            if self.data["asn_pool"] == "none":
                self.data["asn_pool"] = basename(item.meta["pool_file"])
                parsed_name = re.search(
                    _DMS_POOLNAME_REGEX, self.data["asn_pool"].split(".")[0]
                )
                if parsed_name is not None:
                    pool_meta = {
                        "program_id": parsed_name.group(1),
                        "version": parsed_name.group(2),
                    }
                    self.meta["pool_meta"] = pool_meta

        # Product-based updates
        product = self.current_product
        product["name"] = self.dms_product_name

    def make_member(self, item):
        """Create a member from the item

        Parameters
        ----------
        item : dict
            The item to create member from.

        Returns
        -------
        member : Member
            The member
        """
        try:
            exposerr = item["exposerr"]
        except KeyError:
            exposerr = None

        # Get exposure type
        exptype = self.get_exposure_type(item)

        # Determine expected member name
        expname = Utility.rename_to_level2(
            item["filename"], exp_type=item["exp_type"], member_exptype=exptype
        )

        member = Member(
            {
                "expname": expname,
                "exptype": exptype,
                "asn_candidate": item["asn_candidate"],
                "exposerr": exposerr,
            },
            item=item,
        )
        return member

    def make_fov_asn(self):
        """Take the association with an single exposure with _WFI_ in the name
              and expand that to include all 18 detectors.

        Returns
        -------
        associations : [association[, ...]]
            List of new members to be used in place of
            the current one.
        """
        results = []

        # expand the products from _wfi_ to _wfi{det}_
        for product in self["products"]:
            for member in product["members"]:
                asn = copy.deepcopy(self)
                asn.data["products"] = None
                product_name = (
                    splitext(
                        split(self.data["products"][0]["members"][0]["expname"])[1]
                    )[0].rsplit("_", 1)[0]
                    + "_drzl"
                )
                asn.new_product(product_name)
                new_members = asn.current_product["members"]
                if "_wfi_" in member["expname"]:
                    # Make and add a member for each detector
                    for det in WFI_DETECTORS:
                        new_member = copy.deepcopy(member)
                        new_member["expname"] = member["expname"].replace("wfi", det)
                        new_members.append(new_member)
            if asn.is_valid:
                results.append(asn)
                return results
            else:
                return None

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""
        super()._init_hook(item)

        # Set which sequence counter should be used.
        self._sequence = self._sequences[self.data["asn_type"]]

        # Create the product.
        self.new_product()

        # Update meta data
        self.update_asn(item=item)

    def _add(self, item):
        """Add item to this association."""
        member = self.make_member(item)
        if self.is_member(member):
            # logger.debug(
            #     'Member is already part of the association:'
            #     '\n\tassociation: {}'
            #     '\n]tmember: {}'.format(self, member)
            # )
            return

        self.update_validity(member)
        members = self.current_product["members"]
        members.append(member)
        if member["exposerr"] not in _EMPTY:
            logger.warning(
                f"Member {item['filename']} has exposure error \"{member['exposerr']}\""
            )

        # Update meta info
        self.update_asn(item=item, member=member)

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
            self.update_validity(member)
            members.append(member)
        self.sequence = next(self._sequence)

    def __repr__(self):
        # flake8:  noqa: F821
        try:
            file_name, json_repr = self.ioregistry["json"].dump(self)
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


@RegistryMarker.utility
class Utility:
    """Utility functions that understand DMS Level 3 associations"""

    @staticmethod
    def resequence(associations):
        """Resequence the numbering for the ELPP association types"""
        counters = defaultdict(lambda: defaultdict(Counter))
        for asn in associations:
            asn.sequence = next(counters[asn.data["asn_id"]][asn.data["asn_type"]])

    @staticmethod
    def rename_to_level2(level1b_name, exp_type=None, member_exptype="science"):
        """Rename a Level 1b Exposure to a Level2 name.

        The basic transform is changing the suffix `uncal` to
        `cal`, `calints`, or `rate`.

        Parameters
        ----------
        level1b_name : str
            The Level 1b exposure name.

        exp_type:
            JWST exposure type. If not specified,
            it will be presumed that the name
            should get a Level2b name

        is_tso : boolean
            Use 'calints' instead of 'cal' as
            the suffix.

        member_exptype: str
            The association member exposure type, such as "science".

        Returns
        -------
        str
            The Level 2b name
        """
        match = re.match(_LEVEL1B_REGEX, level1b_name)
        if match is None or match.group("type") != "_uncal":
            logger.warning(
                (
                    'Item FILENAME="{}" is not a Level 1b name. Cannot transform to'
                    " Level 2b."
                ).format(level1b_name)
            )
            return level1b_name

        if member_exptype == "background":
            suffix = "x1d"
        else:
            if exp_type in LEVEL2B_EXPTYPES:
                suffix = "cal"
            else:
                suffix = "rate"

        # if is_tso:
        #    suffix += 'ints'

        level2_name = "".join(
            [match.group("path"), "_", suffix, match.group("extension")]
        )
        return level2_name

    @staticmethod
    def get_candidate_list(value):
        """Parse the candidate list from a item value

        Parameters
        ----------
        value : str
            The value from the item to parse. Usually
            item['ASN_CANDIDATE']

        Returns
        -------
        [ACID, ...]
            The list of parsed candidates.
        """

        result = []
        evaled = evaluate(value)
        if is_iterable(evaled):
            result = [ACID(v) for v in evaled]
        return result

    @staticmethod
    @RegistryMarker.callback("finalize")
    def finalize(associations):
        """Check validity and duplications in an association list

        Parameters
        ----------
        associations:[association[, ...]]
            List of associations

        Returns
        -------
        finalized_associations : [association[, ...]]
            The validated list of associations
        """
        finalized_asns = []
        lv3_asns = []
        for asn in associations:
            if isinstance(asn, DMS_ELPP_Base):
                finalized = asn.finalize()
                if finalized is not None:
                    lv3_asns.extend(finalized)
            else:
                finalized_asns.append(asn)

        lv3_asns = prune_duplicate_associations(lv3_asns)
        lv3_asns = prune_duplicate_products(lv3_asns)

        # Ensure sequencing is correct.
        Utility.resequence(lv3_asns)

        # Merge lists and return
        return finalized_asns + lv3_asns


# ---------
# Utilities
# ---------
# Define default product name filling
format_product = FormatTemplate(
    key_formats={"source_id": ["s{:05d}", "s{:s}"], "expspcin": ["{:0>2s}"]}
)


def dms_product_name_noopt(asn):
    """Define product name without any optical elements.

    Parameters
    ---------
    asn : Association
        The association for which the product
        name is to be created.

    Returns
    -------
    product_name : str
        The product name
    """
    target = asn._get_target()

    instrument = asn._get_instrument()

    product_name = "r{}-{}_{}_{}".format(
        #        asn.data['program'],
        asn.data["visi_id"],
        asn.acid.id,
        target,
        instrument,
    )
    return product_name.lower()


def dms_product_name_sources(asn):
    """Produce source-based product names

    Parameters
    ---------
    asn : Association
        The association for which the product
        name is to be created.

    Returns
    -------
    product_name : str
        The product name
    """
    instrument = asn._get_instrument()

    opt_elem = asn._get_opt_element()

    product_name_format = "r{program}-{acid}_{source_id}_{instrument}_{opt_elem}"
    product_name = format_product(
        product_name_format,
        #        program=asn.data['program'],
        program=asn.data["visit_id"],
        acid=asn.acid.id,
        instrument=instrument,
        opt_elem=opt_elem,
    )

    return product_name.lower()


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
class AsnMixin_AuxData:
    """Process special and non-science exposures as science."""

    def get_exposure_type(self, item, default="science"):
        """Override to force exposure type to always be science
        Parameters
        ----------
        item : dict
            The pool entry for which the exposure type is determined
        default : str or None
            The default exposure type.
            If None, routine will raise LookupError
        Returns
        -------
        exposure_type : 'science'
            Returns as science for most Exposures
        exposure_type : 'target_acquisition'
            Returns target_acquisition for mir_tacq
        """
        NEVER_CHANGE = ["target_acquisition"]
        exp_type = super().get_exposure_type(item, default=default)
        if exp_type in NEVER_CHANGE:
            return exp_type
        return "science"


class AsnMixin_Science(DMS_ELPP_Base):
    """Basic science constraints"""

    def __init__(self, *args, **kwargs):
        # Setup target acquisition inclusion
        constraint_acqs = Constraint(
            [
                DMSAttrConstraint(
                    name="acq_obsnum",
                    sources=["obs_num"],
                    value=lambda: "("
                    + "|".join(self.constraints["obs_num"].found_values)
                    + ")",
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

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        self.data["asn_type"] = "spec3"
        super()._init_hook(item)


# ---------------------------------------------
# Mixins to define the broad category of rules.
# ---------------------------------------------
class AsnMixin_Lv2FOV:
    """Level 2 Image association base"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super()._init_hook(item)
        self.data["asn_type"] = "FOV"

    def finalize(self):
        """Finalize association


        Returns
        -------
        associations: [association[, ...]] or None
            List of fully-qualified associations that this association
            represents.
            `None` if a complete association cannot be produced.

        """
        if self.is_valid:
            return self.make_fov_asn()
        else:
            return None


class AsnMixin_Lv2Image:
    """Level 2 Image association base"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super()._init_hook(item)
        self.data["asn_type"] = "image"


class AsnMixin_Lv2GBTDSpass:
    """Level 2 GBTDS association base"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super()._init_hook(item)
        self.data["asn_type"] = "pass"


class AsnMixin_Lv2GBTDSfull:
    """Level 2 GBTDS association base"""

    def _init_hook(self, item):
        """Post-check and pre-add initialization"""

        super()._init_hook(item)
        self.data["asn_type"] = "full"
