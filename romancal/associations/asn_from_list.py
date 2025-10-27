"""Create an association from a list"""

import argparse
import sys
from collections import OrderedDict

from .lib.rules_elpp_base import DMS_ELPP_Base
from .registry import AssociationRegistry

__all__ = ["asn_from_list"]


def asn_from_list(items, rule=DMS_ELPP_Base, **kwargs):
    """Create an association from a list

    Parameters
    ----------
    items: [object [, ...]]
        List of items to add.

    rule: `Association` rule
        The association rule to use.

    kwargs: dict
        Other named parameters required or pertinent to adding
        the items to the association.

    Returns
    -------
    association: `Association`-based instance
        The association with the items added.

    Notes
    -----
    This is a lower-level tool for artificially creating
    an association. As such, the association created may not be valid.
    It is presume the user knows what they are doing.
    """
    asn = rule()
    asn._add_items(items, **kwargs)

    # Preserve specific top-level metadata if provided in kwargs
    if "target" in kwargs.keys():
        target = kwargs["target"]
        asn["target"] = target

    # Always set data_release_id; default to 'p' if not provided
    asn["data_release_id"] = kwargs.get("data_release_id", "p")
    asn = _create_ordered_meta(asn)

    return asn


def _create_ordered_meta(asn):
    """
    Reorder the association metadata so that 'data_release_id' appears
    immediately after 'program' in the top-level dictionary.

    Parameters
    ----------
    asn : Association
        The association object whose metadata should be reordered.

    Returns
    -------
    Association
        The association object with reordered metadata.
    """
    if "program" in asn and "data_release_id" in asn:
        items = list(asn.items())

        # Build new ordered list
        new_items = []
        for k, v in items:
            new_items.append((k, v))
            if k == "program":
                new_items.append(("data_release_id", asn["data_release_id"]))

        # Remove duplicate if present
        seen = set()
        ordered = []
        for k, v in new_items:
            if k == "data_release_id" and k in seen:
                continue
            ordered.append((k, v))
            seen.add(k)

        # Assign to internal dict
        asn.data = OrderedDict(ordered)

    return asn


def _cli(args=None):
    """
    Command-line interface for creating an association from a list of files.

    Parameters
    ----------
    args : list of str or None, optional
        Command line arguments. If None, uses sys.argv[1:].
        If a string, splits on spaces.

    Returns
    -------
    None
    """

    if args is None:
        args = sys.argv[1:]
    if isinstance(args, str):
        args = args.split(" ")

    parser = argparse.ArgumentParser(
        description="Create an association from a list of files",
        usage="asn_from_list -o mosaic_asn.json\n--product-name my_mosaic *.asdf",
    )

    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        required=True,
        help="File to write association to",
    )

    parser.add_argument(
        "--product-name",
        type=str,
        help="The product name when creating a Level 3 association",
    )

    parser.add_argument(
        "--data-release-id",
        type=str,
        default="p",
        help="Data release id to include in the association top-level metadata (default: 'p')",
    )

    parser.add_argument(
        "-r",
        "--rule",
        type=str,
        default="DMS_ELPP_Base",
        help=('The rule to base the association structure on. Default: "%(default)s"'),
    )
    parser.add_argument(
        "--ruledefs",
        action="append",
        help=(
            "Association rules definition file(s) If not specified, the default"
            " rules will be searched."
        ),
    )
    parser.add_argument(
        "-i",
        "--id",
        type=str,
        default="o999",
        help='The association candidate id to use. Default: "%(default)s"',
        dest="acid",
    )
    parser.add_argument(
        "-t",
        "--target",
        type=str,
        default="None",
        help='The target name for the association. Default: "%(default)s"',
        dest="target",
    )

    parser.add_argument(
        "filelist",
        type=str,
        nargs="+",
        help="File list to include in the association",
    )

    parsed = parser.parse_args(args=args)
    print("Parsed args:", parsed)

    # Get the rule
    rule = AssociationRegistry(parsed.ruledefs, include_bases=True)[parsed.rule]

    with open(parsed.output_file, "w") as outfile:
        asn = asn_from_list(
            parsed.filelist,
            rule=rule,
            product_name=parsed.product_name,
            acid=parsed.acid,
            target=parsed.target,
            data_release_id=parsed.data_release_id,
        )
        _, serialized = asn.dump(format="json")
        outfile.write(serialized)
