"""Create an association based on skycells"""

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import roman_datamodels as rdm

import romancal.skycell.match as sm
import romancal.skycell.skymap as sc
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
    output_file_root: os.PathLike | str,
    product_type: str,
    data_release_id: str = "p",
):
    """
    Create the skycell association from a list of L2 calibrated files.

    This function processes a list of Level 2 calibrated files, identifies matching skycells, generates
    associations for skycells. It supports producing associations for three product types:
      - "full": one association per skycell that includes every exposure that touches the skycell.
      - "pass": associations per skycell per pass (subset). The code will loop over passes found in the input.
      - "visit": associations per skycell per visit. The code will loop over visits found in the input.

    Parameters
    ----------
    filelist : list of str
        List of file names to be processed.
    output_file_root : str (or path-like object)
        Root string for the output association file (used by mk_level3_asn_name).
    product_type : str
        Type of product when creating the association (e.g., 'visit', 'pass', 'full').
    data_release_id : str, optional
        Data release identifier to be written into the association top-level metadata.
        Defaults to 'p' when not provided.
    """
    output_file_root = str(output_file_root)

    # normalize product_type
    product_type = (product_type or "full").lower()

    # Split into groups according to requested product_type
    groups = _create_groups(filelist, product_type)
    logger.info(f"Creating {product_type} associations for {len(groups)} groups")

    # For efficiency, precompute intersecting skycells and filter per file
    file_index = _create_intersecting_skycell_index(filelist)

    # Process each group separately
    _process_groups(
        groups,
        file_index,
        output_file_root,
        data_release_id,
        product_type,
    )


def _create_groups(filelist: list[str], product_type: str) -> dict:
    """
    Create groups of files based on the specified product type.

    Parameters
    ----------
    filelist : list of str
        List of filenames to group.
    product_type : str
        Type of product when creating the association (e.g., 'visit', 'pass', 'full').

    Returns
    -------
    dict
        Dictionary mapping group keys to lists of filenames.
    """
    product_type = (product_type or "full").lower()
    if product_type == "visit":
        return _group_files_by_visit(filelist)
    elif product_type == "pass":
        return _group_files_by_pass(filelist)
    else:
        return {"full": list(filelist)}


def _create_intersecting_skycell_index(filelist: list[str]) -> list:
    """
    Create an index of intersecting skycells for each file in the file list.

    Parameters
    ----------
    filelist : list of str
        List of filenames to process.

    Returns
    -------
    list
        List of [filename, skycell_indices, filter_id] records.
    """
    file_index = []
    for file_name in filelist:
        try:
            cal_file = rdm.open(file_name)
            filter_id = cal_file.meta.instrument.optical_element.lower()
            intersecting_skycells = sm.find_skycell_matches(cal_file.meta.wcs)
            cal_file.close()
        except Exception:
            logger.warning(
                "Unable to open %s to read filter or wcs; defaulting to unknown",
                file_name,
            )
            filter_id = "unknown"
            intersecting_skycells = []
        logger.info("Skycell List:%s, %s", file_name, intersecting_skycells)
        file_index.append([file_name, intersecting_skycells, filter_id])
    return file_index


def _process_groups(
    groups: dict,
    file_index: list,
    output_file_root: str,
    data_release_id: str,
    product_type: str,
):
    """
    Process each group of files to create and save skycell associations.

    For each group, this function identifies all unique skycells touched by the group's files,
    determines which files overlap each skycell, constructs the association metadata,
    serializes the association, and writes it to disk.

    Parameters
    ----------
    groups : dict
        Dictionary mapping group keys (e.g., visit or pass identifiers) to lists of filenames.
    file_index : list
        List of [filename, skycell_indices, filter_id] records for all input files.
    output_file_root : str
        Root string for the output association file names.
    data_release_id : str
        Data release identifier to include in the association metadata and filenames.
    product_type : str
        Type of product when creating the association (e.g., 'visit', 'pass', 'full').

    Returns
    -------
    None
    """
    for group_files in groups.values():
        # Build file_list and skycell_indices together
        file_list = [
            [rec[0], rec[1], rec[2]] for rec in file_index if rec[0] in group_files
        ]
        skycell_indices = [idx for rec in file_list for idx in rec[1]]
        # We only want unique skycell indices
        unique_skycell_indices = np.unique(skycell_indices)

        for skycell_index in unique_skycell_indices:
            # Group files by filter for this skycell
            filter_groups = _group_files_by_filter_for_skycell(file_list, skycell_index)
            for filter_id, members in filter_groups.items():
                if not members:
                    continue

                # Get parameters for naming and metadata
                first_member = members[0]
                visit_id_no_r = _extract_visit_id(first_member)
                skycell = sc.SkyCell(skycell_index)
                asn_file_name = mk_level3_asn_name(
                    visit_id_no_r,
                    output_file_root,
                    filter_id,
                    data_release_id,
                    product_type,
                    skycell.name,
                )

                # Create the association metadata
                prompt_product_asn = _create_metadata(
                    members, data_release_id, asn_file_name, skycell, visit_id_no_r
                )

                # Serialize and save the association
                _, serialized = prompt_product_asn.dump(format="json")
                _save_association(asn_file_name, serialized)


def _create_metadata(
    member_list: list[str],
    data_release_id: str,
    asn_file_name: str,
    skycell,
    visit_id_no_r: str,
):
    """
    Create and populate the metadata dictionary for a skycell association.

    Parameters
    ----------
    member_list : list of str
        List of filenames to include in the association.
    data_release_id : str
        Data release identifier to include in the association metadata.
    asn_file_name : str
        Product name for the association, used as the product_name in the ASN.
    skycell : SkyCell
        SkyCell object representing the current skycell.
    visit_id_no_r : str
        Visit ID string (without leading 'r') used for program extraction.

    Returns
    -------
    dict
        Metadata dictionary for the association, ready for serialization.
    """
    prompt_product_asn = asn_from_list.asn_from_list(
        member_list, product_name=asn_file_name
    )
    prompt_product_asn["asn_type"] = "image"
    try:
        program_id = parse_visitID(visit_id_no_r).get("Program", "")
    except Exception:
        program_id = ""
    prompt_product_asn["program"] = program_id
    prompt_product_asn["data_release_id"] = data_release_id
    prompt_product_asn["target"] = skycell.name
    prompt_product_asn["skycell_wcs_info"] = skycell.wcs_info

    return prompt_product_asn


def _save_association(asn_file_name: str, serialized: str):
    """
    Save the association to a file.

    Parameters
    ----------
    asn_file_name : str
        Base name for the association file.
    serialized : str
        Serialized association content.

    Returns
    -------
    None
    """
    out_name = asn_file_name + "_asn.json"
    with open(out_name, "w") as outfile:
        outfile.write(serialized)
    logger.info("Wrote association: %s", out_name)


def _fetch_filter_for(filename: str, file_index) -> str:
    """
    Retrieve the filter ID for a given filename from the precomputed file index.

    Parameters
    ----------
    filename : str
        The filename to look up.
    file_index : list
        List of [filename, skycell_indices, filter_id] records.

    Returns
    -------
    str
        The filter ID associated with the filename, or "unknown" if not found.
    """
    for rec in file_index:
        if rec[0] == filename:
            return rec[2]
    return "unknown"


def _extract_visit_id(filename):
    """
    Extract the visit ID (without leading 'r') from a filename.

    Parameters
    ----------
    filename : str
        The filename to extract the visit ID from.

    Returns
    -------
    str
        The visit ID string without the leading 'r', if present.
    """
    base = Path(filename).stem
    obs = base.split("_")[0]
    if obs.startswith("r"):
        obs = obs[1:]
    return obs[:19]


def _group_files_by_visit(filelist: list[str]) -> dict:
    """
    Group files by Visit_ID (PPPPPCCAAASSSOOOVVV).

    Parameters
    ----------
    filelist : list of str
        List of filenames to group.

    Returns
    -------
    dict
        Dictionary mapping visit_id_no_r (19-char Visit_ID without leading 'r')
        to a list of filenames.
    """
    groups: dict = {}
    for f in filelist:
        visit_id_no_r = _extract_visit_id(f)
        groups.setdefault(visit_id_no_r, []).append(f)
    return groups


def _group_files_by_pass(filelist: list[str]) -> dict:
    """
    Group files by pass identifier (CCAAA portion of Visit_ID).

    Parameters
    ----------
    filelist : list of str
        List of filenames to group.

    Returns
    -------
    dict
        Dictionary mapping pass_key (5 chars CCAAA) to a list of filenames.
    """
    groups: dict = {}
    for f in filelist:
        visit_id_no_r = _extract_visit_id(f)
        # though CCAAA occupies positions 5:10 (0-based indexing),
        # 'program' (PPPPP) is also relevant here, so we include it too
        pass_key = visit_id_no_r[:10] if len(visit_id_no_r) >= 10 else visit_id_no_r
        groups.setdefault(pass_key, []).append(f)
    return groups


def _group_files_by_filter_for_skycell(file_list, skycell_index):
    """
    Group files by filter for a specific skycell index.

    Parameters
    ----------
    file_list : list
        List of [filename, skycell_indices, filter_id] for the group.
    skycell_index : str or int
        The skycell index to filter on.

    Returns
    -------
    dict
        Dictionary mapping filter_id to list of filenames for this skycell.
    """
    filter_groups = {}
    for fname, skycells, filter_id in file_list:
        if np.isin(skycell_index, skycells):
            filter_groups.setdefault(filter_id, []).append(fname)
    return filter_groups


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
        usage="skycell_asn --product-type visit --data-release-id r0 *_cal.asdf -o r512",
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
        help="The product type when creating the association (visit, pass, full)",
    )

    parser.add_argument(
        "--data-release-id",
        type=str,
        default="p",
        help="Data release id to include in the association top-level metadata and filename (default: 'p')",
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
        parsed.data_release_id,
    )
