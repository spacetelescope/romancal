"""Create an association based on skycells"""

import argparse
import logging
import os
import sys

import numpy as np
import roman_datamodels as rdm

import romancal.skycell.match as sm
from romancal.associations import asn_from_list
from romancal.associations.lib.utilities import mk_level3_asn_name
from romancal.lib.basic_utils import parse_visitID as parse_visitID

__all__ = ["skycell_asn"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")


def skycell_asn(
    filelist: list[str],
    output_file_root: os.PathLike,
    product_type: str,
    release_product: str,
):
    """
    Create the skycell association from the list of files.

    This function processes a list of files, identifies matching skycells, generates
    TAN WCS parameters, and creates an association file for the identified skycells.

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
    skycell_indices = []
    file_list = []
    for file_name in filelist:
        cal_file = rdm.open(file_name)
        filter_id = cal_file.meta.instrument.optical_element.lower()
        intersecting_skycells, _ = sm.find_skycell_matches(cal_file.meta.wcs)
        logger.info(f"Skycell List:{file_name}, {intersecting_skycells}")
        file_list.append([file_name, intersecting_skycells])
        skycell_indices.extend(intersecting_skycells)

    unique_skycell_indices = np.unique(skycell_indices)
    for skycell_index in unique_skycell_indices:
        member_list = []
        skycell = sm.SkyCell(skycell_index)
        for a in file_list:
            if np.isin(skycell_index, a[1]):
                member_list.append(a[0])

        # grab all the wcs parameters needed for generate_tan_wcs
        projcell_info = skycell.wcsinfo
        parsed_visit_id = parse_visitID(member_list[0][1:20])
        program_id = parsed_visit_id["Program"]
        asn_file_name = mk_level3_asn_name(
            member_list[0][1:20],
            output_file_root,
            filter_id,
            release_product,
            product_type,
            skycell.name,
        )

        prompt_product_asn = asn_from_list.asn_from_list(
            member_list, product_name=asn_file_name
        )
        prompt_product_asn["asn_type"] = "image"
        prompt_product_asn["program"] = program_id
        prompt_product_asn["target"] = skycell.name
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
