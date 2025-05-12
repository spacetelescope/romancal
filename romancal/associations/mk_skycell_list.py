"""Create a list of skycells for a level 2 file or files based on skycells"""

import argparse
import logging
import os.path
import sys

import numpy as np
import roman_datamodels as rdm

import romancal.skycell.match as sm

__all__ = ["mk_skycell_list"]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")


def mk_skycell_list(output_dir, filelist):
    """
    Create a list of skycell id's based on a list of level 2 files.

    This function processes a list of calibrated detector files,
    identifies matching skycells on the sky and creates
    file with that list of skycells the level 2 file touches.

    Parameters
    ----------
    filelist : list or str
        List of file names to be processed.
    output_dir : str
        Directory for the output skycell file.
    """

    for file_name in filelist:
        input_dir, input_file = os.path.split(file_name)
        cal_file = rdm.open(file_name)
        intersecting_skycells = sm.find_skycell_matches(cal_file.meta.wcs)
        logger.info(f"Skycell List: {file_name}, {intersecting_skycells}")
        output_file_name = os.path.basename(input_file).split(".")[0]
        if not output_dir:
            output_file_name = os.path.join(input_dir, output_file_name)
        else:
            output_file_name = os.path.join(output_dir, output_file_name)

        with open(output_file_name + ".match", "w") as outfile:
            out_string = (
                file_name + " " + np.array2string(intersecting_skycells, separator=",")
            )
            outfile.write(out_string)


def _cli(args=None):
    """Command-line interface for mk_skycell_list

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
        description="Create list of skycells from a level 2 file or files",
        usage="mk_skycell_list *_cal.asdf ",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="",
        help="The optional directory to write the list of skycells",
    )

    parser.add_argument(
        "filelist",
        type=str,
        nargs="+",
        help="Input file list to generate a list of skycells",
    )

    parsed = parser.parse_args(args=args)
    logger.info("Command-line arguments: %s", parsed)
    mk_skycell_list(
        parsed.output_dir,
        parsed.filelist,
    )
