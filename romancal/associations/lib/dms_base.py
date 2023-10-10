"""Association attributes common to DMS-based Rules"""

from romancal.associations.exceptions import AssociationNotValidError
from romancal.associations.lib.acid import ACIDMixin
from romancal.associations.lib.constraint import (
    AttrConstraint,
    Constraint,
    SimpleConstraint,
)
from romancal.associations.lib.counter import Counter
from romancal.associations.lib.utilities import getattr_from_list

__all__ = ["Constraint_TargetAcq", "Constraint_WFSC", "DMSBaseMixin"]

# Default product name
PRODUCT_NAME_DEFAULT = "undefined"

# DMS file name templates
_ASN_NAME_TEMPLATE_STAMP = "r{program}-{acid}_{stamp}_{type}_{sequence:03d}_asn"
_ASN_NAME_TEMPLATE = "r{program}-{acid}_{type}_{sequence:03d}_asn"

# Acquisition and Confirmation images
ACQ_EXP_TYPES = (
    "nrc_taconfirm",
    "nrc_tacq",
)

# Exposure EXP_TYPE to Association EXPTYPE mapping
# flake8: noqa: E241
EXPTYPE_MAP = {
    "nrc_dark": "dark",
    "nrc_flat": "flat",
    "nrc_focus": "engineering",
    "nrc_led": "engineering",
    "nrc_tacq": "target_acquisition",
    "nrc_taconfirm": "target_acquisition",
}

# Coronographic exposures
CORON_EXP_TYPES = ["mir_4qpm", "mir_lyot", "nrc_coron"]

# Roman WFI detectors
WFI_DETECTORS = [
    "wfi01",
    "wfi02",
    "wfi03",
    "wfi04",
    "wfi05",
    "wfi06",
    "wfi07",
    "wfi08",
    "wfi09",
    "wfi10",
    "wfi11",
    "wfi12",
    "wfi13",
    "wfi14",
    "wfi15",
    "wfi16",
    "wfi17",
    "wfi18",
]

# Exposures that get Level2b processing
IMAGE2_SCIENCE_EXP_TYPES = [
    "wfi_image",
]

IMAGE2_NONSCIENCE_EXP_TYPES = [
    "wfi_focus",
]
IMAGE2_NONSCIENCE_EXP_TYPES.extend(ACQ_EXP_TYPES)

SPEC2_SCIENCE_EXP_TYPES = [
    "wfi_grism",
    "wfi_prism",
]

SPECIAL_EXPOSURE_MODIFIERS = {
    "background": ["bkgdtarg"],
    "imprint": ["is_imprt"],
    "psf": ["is_psf"],
}

#
# Key that uniquely identifies members.
MEMBER_KEY = "expname"

# Non-specified values found in DMS Association Pools
_EMPTY = (
    None,
    "",
    "NULL",
    "Null",
    "null",
    "--",
    "N",
    "n",
    "F",
    "f",
    "FALSE",
    "false",
    "False",
    "N/A",
    "n/a",
)

# Degraded status information
_DEGRADED_STATUS_OK = "No known degraded exposures in association."
_DEGRADED_STATUS_NOTOK = (
    "One or more members have an error associated with them.\nDetails can be found in"
    " the member.exposerr attribute."
)


class DMSBaseMixin(ACIDMixin):
    """Association attributes common to DMS-based Rules

    Attributes
    ----------
    sequence : int
        The sequence number of the current association
    """

    # Associations of the same type are sequenced.
    _sequence = Counter(start=1)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._acid = None
        self._asn_name = None
        self.sequence = None
        if "degraded_status" not in self.data:
            self.data["degraded_status"] = _DEGRADED_STATUS_OK
        if "program" not in self.data:
            self.data["program"] = "noprogram"

    @classmethod
    def create(cls, item, version_id=None):
        """Create association if item belongs

        Parameters
        ----------
        item : dict
            The item to initialize the association with.

        version_id : str or None
            Version_Id to use in the name of this association.
            If None, nothing is added.

        Returns
        -------
        (association, reprocess_list)
            2-tuple consisting of:

                - association : The association or, if the item does not
                  match this rule, None
                - [ProcessList[, ...]]: List of items to process again.
        """
        asn, reprocess = super().create(item, version_id)
        if not asn:
            return None, reprocess
        asn.sequence = next(asn._sequence)
        return asn, reprocess

    @property
    def acid(self):
        """Association ID"""
        acid = self._acid
        if self._acid is None:
            acid = self.acid_from_constraints()
        return acid

    @property
    def asn_name(self):
        """The association name

        The name that identifies this association. When dumped,
        will form the basis for the suggested file name.

        Typically, it is generated based on the current state of
        the association, but can be overridden.
        """
        if self._asn_name:
            return self._asn_name

        program = self.data["program"]
        version_id = self.version_id
        asn_type = self.data["asn_type"]
        sequence = self.sequence

        if version_id:
            name = _ASN_NAME_TEMPLATE_STAMP.format(
                program=program,
                acid=self.acid.id,
                stamp=version_id,
                type=asn_type,
                sequence=sequence,
            )
        else:
            name = _ASN_NAME_TEMPLATE.format(
                program=program,
                acid=self.acid.id,
                type=asn_type,
                sequence=sequence,
            )
        return name.lower()

    @asn_name.setter
    def asn_name(self, name):
        """Override calculated association name"""
        self._asn_name = name

    @property
    def current_product(self):
        """Return currnet products"""
        return self.data["products"][-1]

    @property
    def from_items(self):
        """The list of items that contributed to the association."""
        try:
            items = [
                member.item
                for product in self["products"]
                for member in product["members"]
            ]
        except KeyError:
            items = []
        return items

    @property
    def member_ids(self):
        """Set of all member ids in all products of this association"""
        member_ids = {
            member[MEMBER_KEY]
            for product in self["products"]
            for member in product["members"]
        }
        return member_ids

    @property
    def validity(self):
        """Keeper of the validity tests"""
        try:
            validity = self._validity
        except AttributeError:
            self._validity = {}
            validity = self._validity
        return validity

    @validity.setter
    def validity(self, item):
        """Set validity dict"""
        self._validity = item

    def get_exposure_type(self, item, default="science"):
        """Determine the exposure type of a pool item

        Parameters
        ----------
        item : dict
            The pool entry to determine the exposure type of

        default : str or None
            The default exposure type.
            If None, routine will raise LookupError

        Returns
        -------
        exposure_type : str
            Exposure type. Can be one of

                - 'science': Item contains science data
                - 'target_acquisition': Item contains target acquisition data.
                - 'autoflat': NIRSpec AUTOFLAT
                - 'autowave': NIRSpec AUTOWAVE
                - 'psf': PSF
                - 'imprint': MSA/IFU Imprint/Leakcal

        Raises
        ------
        LookupError
            When `default` is None and an exposure type cannot be determined
        """
        return get_exposure_type(item, default=default, association=self)

    def is_member(self, new_member):
        """Check if member is already a member

        Parameters
        ----------
        new_member : Member
            The member to check for
        """
        try:
            current_members = self.current_product["members"]
        except KeyError:
            return False

        for member in current_members:
            if member == new_member:
                return True
        return False

    def is_item_member(self, item):
        """Check if item is already a member of this association

        Parameters
        ----------
        item : dict
            The item to check for.

        Returns
        -------
        is_item_member : bool
            True if item is a member.
        """
        return item in self.from_items

    def item_getattr(self, item, attributes):
        """Return value from any of a list of attributes

        Parameters
        ----------
        item : dict
            item to retrieve from

        attributes : list
            List of attributes

        Returns
        -------
        (attribute, value)
            Returns the value and the attribute from
            which the value was taken.

        Raises
        ------
        KeyError
            None of the attributes are found in the dict.
        """
        return item_getattr(item, attributes, self)

    def new_product(self, product_name=PRODUCT_NAME_DEFAULT):
        """Start a new product"""
        product = {"name": product_name, "members": []}
        try:
            self.data["products"].append(product)
        except (AttributeError, KeyError):
            self.data["products"] = [product]

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
        self.update_degraded_status()
        self.set_visit_id(item)

    def set_visit_id(self, item):
        """Set the visit id in the association"""

        if "visit_id" in self.data:
            pass
        else:
            self.data["visit_id"] = item["visit_id"]
            try:
                self.data["visit_id"] = item["visit_id"]
            except KeyError:
                # logger.debug(f'Visit_id not found')
                self.data["visit_id"] = "None"

    def update_degraded_status(self):
        """Update association degraded status"""

        if self.data["degraded_status"] == _DEGRADED_STATUS_OK:
            for product in self.data["products"]:
                for member in product["members"]:
                    try:
                        exposerr = member["exposerr"]
                    except KeyError:
                        continue
                    else:
                        if exposerr not in _EMPTY:
                            self.data["degraded_status"] = _DEGRADED_STATUS_NOTOK
                            break

    def update_validity(self, entry):
        for test in self.validity.values():
            if not test["validated"]:
                test["validated"] = test["check"](entry)

    @classmethod
    def reset_sequence(cls):
        cls._sequence = Counter(start=1)

    @classmethod
    def validate(cls, asn):
        super().validate(asn)

        if isinstance(asn, DMSBaseMixin):
            result = False
            try:
                result = all(test["validated"] for test in asn.validity.values())
            except (AttributeError, KeyError):
                raise AssociationNotValidError("Validation failed")
            if not result:
                raise AssociationNotValidError("Validation failed validity tests.")
        return True

    def _get_exposure(self):
        """Get string representation of the exposure id

        Returns
        -------
        exposure : str
            The Level3 Product name representation
            of the exposure & activity id.
        """
        exposure = ""
        try:
            activity_id = format_list(self.constraints["activity_id"].found_values)
        except KeyError:
            pass
        else:
            if activity_id not in _EMPTY:
                exposure = f"{activity_id:0>2s}"
        return exposure

    def _get_instrument(self):
        """Get string representation of the instrument

        Returns
        -------
        instrument : str
            The Level3 Product name representation
            of the instrument
        """
        instrument = format_list(self.constraints["instrument"].found_values)
        return instrument

    def _get_opt_element(self):
        """Get string representation of the optical elements

        Returns
        -------
        opt_elem : str
            The Level3 Product name representation
            of the optical elements.
        """
        # Retrieve all the optical elements
        opt_elems = []
        for opt_elem in ["opt_elem", "opt_elem2", "opt_elem3"]:
            try:
                values = list(self.constraints[opt_elem].found_values)
            except KeyError:
                pass
            else:
                values.sort(key=str.lower)
                value = format_list(values)
                if value not in _EMPTY:
                    opt_elems.append(value)

        # Build the string. Sort the elements in order to
        # create data-independent results
        opt_elems.sort(key=str.lower)
        opt_elem = "-".join(opt_elems)
        if opt_elem == "":
            opt_elem = "clear"

        return opt_elem

    def _get_subarray(self):
        """Get string representation of the subarray

        Returns
        -------
        subarray : str
            The Level3 Product name representation
            of the subarray.
        """
        result = ""
        try:
            subarray = format_list(self.constraints["seq_id"].found_values)
        except KeyError:
            subarray = None
            return result
        if subarray == 0:
            subarray = None
        try:
            format_list(self.constraints["subcat"].found_values)
        except KeyError:
            pass

        return result

    def _get_target(self):
        """Get string representation of the target

        Returns
        -------
        target : str
            The Level3 Product name representation
            of the target or source ID.
        """
        target_id = format_list(self.constraints["target"].found_values)
        target = f"t{str(target_id):0>3s}"
        return target

    def _get_grating(self):
        """Get string representation of the grating in use

        Returns
        -------
        grating : str
            The Level3 Product name representation
            of the grating in use.
        """
        grating_id = format_list(self.constraints["grating"].found_values)
        grating = f"{str(grating_id):0>3s}"
        return grating


# -----------------
# Basic constraints
# -----------------
class DMSAttrConstraint(AttrConstraint):
    """DMS-focused attribute constraint

    Forces definition of invalid values
    """

    def __init__(self, **kwargs):
        if kwargs.get("invalid_values", None) is None:
            kwargs["invalid_values"] = _EMPTY

        super().__init__(**kwargs)


class Constraint_TargetAcq(SimpleConstraint):
    """Select on target acquisition exposures

    Parameters
    ----------
    association:  ~romancal.associations.Association
        If specified, use the `get_exposure_type` method
        of the association rather than the utility version.
    """

    def __init__(self, association=None):
        if association is None:
            _get_exposure_type = get_exposure_type
        else:
            _get_exposure_type = association.get_exposure_type

        super().__init__(
            name="target_acq", value="target_acquisition", sources=_get_exposure_type
        )


class Constraint_WFSC(Constraint):
    """Match on Wave Front Sensing and Control Observations"""

    def __init__(self, *args, **kwargs):
        super().__init__(
            [
                Constraint(
                    [
                        DMSAttrConstraint(
                            name="wfsc",
                            sources=["visitype"],
                            value=".+wfsc.+",
                            force_unique=True,
                        )
                    ]
                )
            ]
        )


# #########
# Utilities
# #########
def format_list(alist):
    """Format a list according to DMS naming specs"""
    return "-".join(alist)


def get_exposure_type(item, default="science", association=None):
    """Determine the exposure type of a pool item

    Parameters
    ----------
    item : dict
        The pool entry to determine the exposure type of

    default : str or None
        The default exposure type.
        If None, routine will raise LookupError



    Returns
    -------
    exposure_type : str
        Exposure type. Can be one of

        - 'science': Item contains science data
        - 'target_acquisition': Item contains target acquisition data.
        - 'autoflat': NIRSpec AUTOFLAT
        - 'autowave': NIRSpec AUTOWAVE
        - 'psf': PSF
        - 'imprint': MSA/IFU Imprint/Leakcal

    Raises
    ------
    LookupError
        When `default` is None and an exposure type cannot be determined
    """

    # Specify how attributes of the item are retrieved.
    def _item_attr(item, sources):
        """Get attribute value of an item

        This simplifies the call to `item_getattr`
        """
        source, value = item_getattr(item, sources, association=association)
        return value

    # Define default type.
    result = default

    # Retrieve pointing type. This decides the basic exposure type.
    # If the pointing is not science, we're done.
    try:
        result = _item_attr(item, ["pntgtype"])
    except KeyError:
        pass
    else:
        if result != "science":
            return result

    # We have a science exposure. Refine further.
    #
    # Base type off of exposure type.
    try:
        exp_type = _item_attr(item, ["exp_type"])
    except KeyError:
        raise LookupError("Exposure type cannot be determined")

    result = EXPTYPE_MAP.get(exp_type, default)

    if result is None:
        raise LookupError("Cannot determine exposure type")

    # If result is not science, we're done.
    if result != "science":
        return result

    # For `science` data, compare against special modifiers
    # to further refine the type.
    for special, source in SPECIAL_EXPOSURE_MODIFIERS.items():
        try:
            _item_attr(item, source)
        except KeyError:
            pass
        else:
            result = special
            break

    return result


def item_getattr(item, attributes, association=None):
    """Return value from any of a list of attributes

    Parameters
    ----------
    item : dict
        item to retrieve from

    attributes : list
        List of attributes

    Returns
    -------
    (attribute, value)
        Returns the value and the attribute from
        which the value was taken.

    Raises
    ------
    KeyError
        None of the attributes are found in the dict.
    """
    if association is None:
        invalid_values = _EMPTY
    else:
        invalid_values = association.INVALID_VALUES
    return getattr_from_list(item, attributes, invalid_values=invalid_values)
