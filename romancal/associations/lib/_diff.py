""" "Diff and compare associations"""

import logging
from collections import UserList
from copy import copy

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# #########################
# Define the types of diffs
# #########################
class DiffError(AssertionError):
    """Base Class for difference errors"""


class MemberMismatchError(DiffError):
    """Membership does not match"""


class MultiDiffError(UserList, DiffError):
    """List of diff errors"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self):
        message = ["Following diffs found:\n"]
        for diff in self:
            message.extend(["\n****\n", str(diff), "\n"])
        return "".join(message)


def compare_product_membership(left, right):
    """Compare membership between products

    Parameters
    ---------
    left, right : dict
        Two, individual, association products to compare

    Raises
    ------
    MultiDiffError
        If there are differences. The message will contain
        all the differences.
    """
    diffs = MultiDiffError()
    if len(right["members"]) != len(left["members"]):
        diffs.append(
            MemberMismatchError(
                "Product Member length differs:"
                " Left Product {left_product_name} len {left_len} !=  "
                " Right Product {right_product_name} len {right_len}"
                "".format(
                    left_product_name=left["name"],
                    left_len=len(left["members"]),
                    right_product_name=right["name"],
                    right_len=len(right["members"]),
                )
            )
        )

    members_right = copy(right["members"])
    for left_member in left["members"]:
        for right_member in members_right:
            if left_member["expname"] != right_member["expname"]:
                continue

            if left_member["exptype"] != right_member["exptype"]:
                diffs.append(
                    MemberMismatchError(
                        "Left {left_expname}:{left_exptype} != Right"
                        " {right_expname}:{right_exptype}".format(
                            left_expname=left_member["expname"],
                            left_exptype=left_member["exptype"],
                            right_expname=right_member["expname"],
                            right_exptype=right_member["exptype"],
                        )
                    )
                )

            members_right.remove(right_member)
            break
        else:
            diffs.append(
                MemberMismatchError(
                    f"Left {left_member['expname']}:{left_member['exptype']} has no"
                    " counterpart in right"
                )
            )

    if len(members_right) != 0:
        diffs.append(
            MemberMismatchError(
                "Right has {len_over} unaccounted for members starting with"
                " {right_expname}:{right_exptype}".format(
                    len_over=len(members_right),
                    right_expname=members_right[0]["expname"],
                    right_exptype=members_right[0]["exptype"],
                )
            )
        )

    if diffs:
        raise diffs
