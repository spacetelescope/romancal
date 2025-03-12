"""Create an association based on skycells"""

import argparse
import logging
import sys

import numpy as np
import roman_datamodels as rdm

import romancal.patch_match.patch_match as pm
from romancal.patch_match.patch_match import get_projectioncell_wcs
from romancal.associations import asn_from_list
from romancal.lib.basic_utils import parse_visitID as parse_visitID

__all__ = ["skycell_asn"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")


def skycell_asn(filelist, output_file_root, product_type, release_product):
    """
    Create the skycell association from the list of files.

    This function processes a list of files, identifies matching patches, generates
    TAN WCS parameters, and creates an association file for the identified sky patches.

    Parameters
    ----------
    filelist : list of str
        List of file names to be processed.
    output_file_root : str
        Root string for the output association file.
    product_type : str
        Type of product when creating the association (e.g., 'visit', 'daily').
    release_product : str
        Release product identifier for the association.

    Returns
    -------
    None
    """
    all_patches = []
    file_list = []
    for file_name in filelist:
        cal_file = rdm.open(file_name)
        filter_id = cal_file.meta.instrument.optical_element.lower()
        file_patch_list = pm.find_patch_matches(cal_file.meta.wcs)
        logger.info(f"Patch List:{file_name}, {file_patch_list[0]}")
        file_list.append([file_name, file_patch_list[0]])
        all_patches.append(file_patch_list[0])

    unique_patches = np.unique(np.concatenate(all_patches))
    for item in unique_patches:
        member_list = []
        patch_name = pm.PATCH_TABLE[item]["name"]
        for a in file_list:
            if np.isin(item, a[1]):
                member_list.append(a[0])

        # grab all the wcs parameters needed for generate_tan_wcs
        projcell_info = get_projectioncell_wcs(item)
        parsed_visit_id = parse_visitID(member_list[0][1:20])
        program_id = parsed_visit_id["Program"]
        root_asn_name = output_file_root
        product_release = release_product
        sep = "_"

        product_name_mapping = {
            "visit": "v"
            + parsed_visit_id["Execution"]
            + parsed_visit_id["Pass"]
            + parsed_visit_id["Segment"]
            + parsed_visit_id["Observation"],
            "daily": "d"
            + parsed_visit_id["Execution"]
            + parsed_visit_id["Pass"]
            + parsed_visit_id["Segment"],
            "pass": "p" + parsed_visit_id["Execution"] + parsed_visit_id["Pass"],
            "full": "full",
            "user": "user",
        }

        pr_name = product_name_mapping.get(product_type, "unknown")

        asn_file_name = (
            root_asn_name
            + sep
            + product_release
            + sep
            + pr_name
            + sep
            + patch_name
            + sep
            + filter_id
        )

        prompt_product_asn = asn_from_list.asn_from_list(
            member_list, product_name=asn_file_name
        )
        prompt_product_asn["asn_type"] = "image"
        prompt_product_asn["program"] = program_id
        prompt_product_asn["target"] = patch_name
        prompt_product_asn["skycell_wcs_info"] = projcell_info

        _, serialized = prompt_product_asn.dump(format="json")

        with open(asn_file_name + "_asn.json", "w") as outfile:
            outfile.write(serialized)


def _cli(args=None):
    """Command-line interface for list_to_asn

    Parameters
    ----------
    args: [str, ...], or None
        The command line arguments. Can be one of
            - `None`: `sys.argv` is then used.
            - `[str, ...]`: A list of strings which create the command line
              with the similar structure as `sys.argv`
    """
    if args is None:
        args = sys.argv[1:]
    if isinstance(args, str):
        args = args.split(" ")

    parser = argparse.ArgumentParser(
        description="Create an association from a list of files",
        usage="skycell_asn --product-type visit --release-product prompt *_cal.asdf -o r512",
    )

    parser.add_argument(
        "-o",
        "--output-file-root",
        type=str,
        required=True,
        help="Root string for file to write association to",
    )

    parser.add_argument(
        "-f",
        "--format",
        type=str,
        default="json",
        help='Format of the association files. Default: "%(default)s"',
    )

    parser.add_argument(
        "--product-type",
        type=str,
        default="visit",
        help="The product type when creating the association",
    )

    parser.add_argument(
        "--release-product",
        type=str,
        default="p",
        help="The release product when creating the association",
    )

    parser.add_argument(
        "-r",
        "--rule",
        type=str,
        default="DMS_ELPP_Base",
        help=('The rule to base the association structure on. Default: "%(default)s"'),
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
        "filelist",
        type=str,
        nargs="+",
        help="File list to include in the association",
    )

    parsed = parser.parse_args(args=args)
    logger.info("Command-line arguments: %s", parsed)

    skycell_asn(
        parsed.filelist,
        parsed.output_file_root,
        parsed.product_type,
        parsed.release_product,
    )
