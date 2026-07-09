"""Association attributes common to DMS-based Rules"""

__all__ = ["DMSBaseMixin"]

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


class DMSBaseMixin:
    """Association attributes common to DMS-based Rules"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._asn_name = None
        if "degraded_status" not in self.data:
            self.data["degraded_status"] = _DEGRADED_STATUS_OK
        if "program" not in self.data:
            self.data["program"] = "noprogram"

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
        # sequence was a class attribute incremented several times based on test order
        sequence = 1
        # acidid was always a3001
        acidid = "a3001"
        target = self.target

        if version_id:
            name = _ASN_NAME_TEMPLATE_STAMP.format(
                program=program,
                acid=acidid,
                stamp=version_id,
                type=asn_type,
                sequence=sequence,
                target=target,
            )
        else:
            name = _ASN_NAME_TEMPLATE.format(
                program=program,
                acid=acidid,
                type=asn_type,
                sequence=sequence,
                target=target,
            )
        return name.lower()

    @asn_name.setter
    def asn_name(self, name):
        """Override calculated association name"""
        self._asn_name = name

    @property
    def member_ids(self):
        """Set of all member ids in all products of this association"""
        member_ids = {
            member[MEMBER_KEY]
            for product in self["products"]
            for member in product["members"]
        }
        return member_ids

    def new_product(self, product_name=PRODUCT_NAME_DEFAULT):
        """Start a new product"""
        product = {"name": product_name, "members": []}
        try:
            self.data["products"].append(product)
        except (AttributeError, KeyError):
            self.data["products"] = [product]
