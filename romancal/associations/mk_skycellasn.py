"""Create a set of level 3 associations based on skycells"""

import argparse
import logging
import sys
from os.path import basename
import pdb

import numpy as np

import romancal.patch_match.patch_match as pm
from romancal.associations import asn_from_list
from romancal.lib.basic_utils import parse_visitID as parse_visitID


__all__ = ["mk_skycellasn"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")


def mk_skycellasn(filelist, release_product, optical_element):
    """
    Create level 3 associations from a list of match files generated with mk_patchlist.

    This function processes a list of match files, sorts them based on which input exposures
    contribute to a given skycell and output an association with the needed exposure files. 

    Parameters
    ----------
    filelist : list or str
        List of file names to be processed.
    """

    file_list = []
    patch_list = []
    for file_name in filelist:
        with open(file_name) as match_file:
            logger.debug(f"Working on file:{file_name}")
            match_list = match_file.read().split(" ")
            match_list[1]= np.fromstring(match_list[1].strip('[]'), dtype='int', sep=',')
            patch_list.append(match_list)


    # get a list of all patches in the match files
    patch_array = []
    for item in patch_list:
        patch_array.append(item[1])
        #pdb.set_trace()

    logger.debug(f"patch_array: {patch_array}")
    unique_patches = np.unique(np.concatenate(patch_array))
    logger.info(f"Unique Patches: {unique_patches}")
    pm.load_patch_table()
    for item in unique_patches:
        member_list = []
        patch_name = pm.PATCH_TABLE[item]["name"]
        for entry in patch_list:
            if np.isin(item, entry[1]):
                member_list.append(entry[0])

        # grab all the wcs parameters needed for generate_tan_wcs
        projcell_info = get_projectioncell_wcs(item)
        parsed_visit_id = parse_visitID(member_list[0][1:20])
        program_id = parsed_visit_id["Program"]

        asn_file_name = ('r' + program_id + '_'
                         + release_product + '_'
                         + 'v' + member_list[0][1:20] + '_'
                         + patch_name +  '_'
                         + optical_element)
        prompt_product_asn = asn_from_list.asn_from_list(
            member_list, product_name=asn_file_name
        )
        prompt_product_asn["asn_type"] = "image"
        prompt_product_asn["program"] = program_id
        prompt_product_asn["target"] = patch_name
        prompt_product_asn["skycell_wcs_info"] = projcell_info

        _, serialized = prompt_product_asn.dump(format="json")

        logger.info(f"Writing association with root name: {asn_file_name}")

        with open(asn_file_name + "_asn.json", "w") as outfile:
            outfile.write(serialized)

def get_projectioncell_wcs(index):
    """Return the projection cell wcs info as a dictionary based on the db index number"""

    # check to see if an index is being passed
    if not isinstance(index, np.int64):
        logger.info(f"Input index needs to be a numpy int64 variable")
        return None
    
    # check to see if the patch table is loaded if not load it
    if pm.PATCH_TABLE is None:
        pm.load_patch_table()
        
    projcell_info = dict(
        [
            ("name", pm.PATCH_TABLE[index]["name"]),
            ("pixel_scale", float(pm.PATCH_TABLE[index]["pixel_scale"])),
            (
                "ra_projection_center",
                float(pm.PATCH_TABLE[index]["ra_projection_center"]),
            ),
            (
                "dec_projection_center",
                float(pm.PATCH_TABLE[index]["dec_projection_center"]),
            ),
            ("x0_projection", float(pm.PATCH_TABLE[index]["x0_projection"])),
            ("y0_projection", float(pm.PATCH_TABLE[index]["y0_projection"])),
            ("ra_center", float(pm.PATCH_TABLE[index]["ra_center"])),
            ("dec_center", float(pm.PATCH_TABLE[index]["dec_center"])),
            ("nx", int(pm.PATCH_TABLE[index]["nx"])),
            ("ny", int(pm.PATCH_TABLE[index]["ny"])),
            ("orientat", float(pm.PATCH_TABLE[index]["orientat"])),
            (
                "orientat_projection_center",
                float(pm.PATCH_TABLE[index]["orientat_projection_center"]),
            ),
        ]
    )

    return projcell_info

def _cli(args=None):
    """Command-line interface for mk_skycellasn

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
        usage="mk_skycellasn *.match ",
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
        "--optical_element",
        type=str,
        default="f158",
        help="The optical element used for the visit",
        dest="optical_element"
    )

    parsed = parser.parse_args(args=args)
    logger.info("Command-line arguments: %s", parsed)
    mk_skycellasn(
        parsed.filelist,
        parsed.release_product,
        parsed.optical_element,
    )
    
