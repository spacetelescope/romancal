#!/usr/bin/env python

"""Add velocity aberration correction information to the FITS files provided."""

import argparse
import sys

from roman_datamodels import open as dm_open
from stcal.velocity_aberration import compute_va_effects

__all__ = []  # type: ignore[var-annotated]


def parse_args(args):
    """
    Parse given arguments.

    Parameters
    ----------
    args : command-line input
        The arguments to parse

    Returns
    -------
    parser : object
        Parsed arguments
    """
    description_text = """
    Add velocity aberration correction information to Roman data.
    Input via command line can be one or more files.
    It assumes the following meta keys are present:
    meta.ephemeris.velocity_x (km/sec),
    meta.ephemeris.velocity_y (km/sec),
    meta.ephemeris.velocity_z (km/sec),
    meta.wcsinfo.dec_ref (deg),
    meta.wcsinfo.ra_ref (deg),

    The meta keys updates are:
    meta.velocity_aberration.dec_reference (deg)
    meta.velocity_aberration.ra_reference (deg)
    meta.velocity_aberration.scale_factor (unitless multiplication factor)
    """
    parser = argparse.ArgumentParser(
        prog="set_velocity_aberration",
        description=description_text,
    )
    parser.add_argument(
        "filename",
        nargs="+",
        help="Roman file(s) to which to add the velocity aberration correction information.",
    )
    return parser.parse_args(args)


def main():
    """Parse arguments and add velocity aberration correction information to the files provided."""
    args = parse_args(sys.argv[1:])
    for filename in args.filename:
        add_dva(filename)


def add_dva(filename):
    """
    Determine velocity aberration.

    Given the name of a valid partially populated level ScienceRawModel file,
    determine the velocity aberration scale factor and apparent target position
    in the moving (telescope) frame.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    filename : str
        The name of the file to be updated.
    """
    model = dm_open(filename)
    scale_factor, apparent_ra, apparent_dec = compute_va_effects(
        velocity_x=model.meta.ephemeris.velocity_x,
        velocity_y=model.meta.ephemeris.velocity_y,
        velocity_z=model.meta.ephemeris.velocity_z,
        ra=model.meta.wcsinfo.ra_ref,
        dec=model.meta.wcsinfo.dec_ref,
    )

    # update header
    model.meta.velocity_aberration.scale_factor = scale_factor
    model.meta.velocity_aberration.ra_reference = apparent_ra
    model.meta.velocity_aberration.dec_reference = apparent_dec
    model.save(filename)


if __name__ == "__main__":
    main()
