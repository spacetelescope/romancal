"""Create a list of patches for a level 2 file or files based on skycells"""

import argparse
import logging
import sys
from os.path import basename

import numpy as np
import roman_datamodels as rdm

import romancal.patch_match.patch_match as pm
from romancal.lib.basic_utils import parse_visitID as parse_visitID

__all__ = ["mk_patchlist"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")


def mk_patchlist(filelist):
    """
    Create a list of skycell id's based on a list of level 2 files.

    This function processes a list of calibrated detector files,
    identifies matching patches on the sky and creates
    file with that list of sky patches the level 2 file touches.

    Parameters
    ----------
    filelist : list or str
        List of file names to be processed.
    """

    all_patches = []
    file_list = []
    for file_name in filelist:
        cal_file = rdm.open(file_name)
        file_patch_list = pm.find_patch_matches(cal_file.meta.wcs)
        logger.info(f"Patch List:{file_name}, {file_patch_list[0]}")
        output_file_name = basename(file_name).split(".")[0]
        with open(output_file_name + ".match", "w") as outfile:
            out_string = file_name+ ' ' + np.array2string(file_patch_list[0], separator=',')
            outfile.write(out_string)
        

def _cli(args=None):
    """Command-line interface for mk_patchlist

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
        description="Create list of patches from a level 2 file or files",
        usage="mk_patchlist *_cal.asdf ",
    )
    parser.add_argument(
        "filelist",
        type=str,
        nargs="+",
        help="Input file list to generate a list of patchs",
    )

    parsed = parser.parse_args(args=args)
    logger.info("Command-line arguments: %s", parsed)
    mk_patchlist(
        parsed.filelist,
    )
