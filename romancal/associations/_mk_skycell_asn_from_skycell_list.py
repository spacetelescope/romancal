"""Create a set of level 3 associations based on skycells"""

import argparse
import logging
import os.path
import sys

import numpy as np

import romancal.skycell.skymap as sc
from romancal.associations import asn_from_list
from romancal.associations.lib.utilities import mk_level3_asn_name
from romancal.lib.basic_utils import parse_visitID

__all__ = ["mk_skycell_asn_from_skycell_list"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")


def mk_skycell_asn_from_skycell_list(
    filelist, release_product, product_type, optical_element
):
    """
    Create level 3 associations from a list of match files generated with mk_skycell_list.

    This function processes a list of match files, sorts them based on which input exposures
    contribute to a given skycell and output an association with the needed exposure files.

    Parameters
    ----------
    filelist : list or str
        List of file names to be processed.
    release-product: str
        The name of the release product
    optical_element: str
        The name of the optical element for the associations
    """

    skycell_indices = []
    for file_name in filelist:
        with open(file_name) as match_file:
            logger.debug(f"Working on file:{file_name}")
            match_list = match_file.read().split(" ")
            match_list[1] = np.fromstring(
                match_list[1].strip("[]"), dtype="int", sep=","
            )
            skycell_indices.append(match_list)

    # get a list of all skycells in the match files
    skycell_indices_array = [item[1] for item in skycell_indices]

    logger.debug(f"All Skycells: {skycell_indices_array}")
    skycells = sc.SkyCells(np.unique(np.concatenate(skycell_indices_array)))
    logger.info(f"Unique Skycells: {skycells}")
    for skycell_index, skycell_name, skycell_wcs_info in zip(
        skycells.indices, skycells.names, skycells.wcs_infos, strict=True
    ):
        member_list = []
        for entry in skycell_indices:
            if np.isin(skycell_index, entry[1]):
                member_list.append(os.path.basename(entry[0]).split(".")[0])

        # grab all the wcs parameters needed for generate_tan_wcs
        parsed_visit_id = parse_visitID(member_list[0][1:20])
        program_id = parsed_visit_id["Program"]
        output_file_root = "r" + program_id
        asn_file_name = mk_level3_asn_name(
            member_list[0][1:20],
            output_file_root,
            optical_element,
            release_product,
            product_type,
            skycell_name,
        )

        prompt_product_asn = asn_from_list.asn_from_list(
            member_list, product_name=asn_file_name
        )
        prompt_product_asn["asn_type"] = "image"
        prompt_product_asn["program"] = program_id
        prompt_product_asn["target"] = skycell_name
        prompt_product_asn["skycell_wcs_info"] = skycell_wcs_info

        _, serialized = prompt_product_asn.dump(format="json")

        logger.info(f"Writing association with root name: {asn_file_name}")

        with open(asn_file_name + "_asn.json", "w") as outfile:
            outfile.write(serialized)


def _cli(args=None):
    """Command-line interface for mk_skycell_asn_from_list

    Parameters
    ----------
    args: [str, ...], or None
        The command line arguments. Can be one of
            - `None`: `sys.argv` is then used.
            - `[str, ...]`: A list of strings which create the command line
              with the similar structure as `sys.argv`
    """

    def __init__(self, args=None):
        self.configure(args)

    if args is None:
        args = sys.argv[1:]
    if isinstance(args, str):
        args = args.split(" ")

    parser = argparse.ArgumentParser(
        description="Create level 3 associations from a list of match files",
        usage="mk_skycell_asn_from_skycell_list *.match ",
    )
    parser.add_argument(
        "filelist",
        type=str,
        nargs="+",
        help="A list of match files to generate level 3 asn's",
    )

    parser.add_argument(
        "--release-product",
        type=str,
        default="p",
        help="The release product when creating the association",
    )
    parser.add_argument(
        "--product-type",
        type=str,
        default="visit",
        help="The product type when creating the association (visit, pass, ....",
    )
    parser.add_argument(
        "--optical_element",
        type=str,
        default="f158",
        help="The optical element used for the visit",
        dest="optical_element",
    )

    parsed = parser.parse_args(args=args)
    logger.info("Command-line arguments: %s", parsed)
    mk_skycell_asn_from_skycell_list(
        parsed.filelist,
        parsed.release_product,
        parsed.product_type,
        parsed.optical_element,
    )
