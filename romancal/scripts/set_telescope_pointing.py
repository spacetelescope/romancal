#!/usr/bin/env python

"""
Set the initial world coordinate system for Roman exposures.

The Roman engineering database is queried for the Roman Space Telescope observatory
orientation parameters, and converts that orientation to a WCS
for a list of exposures.
"""

# Licensed under a 3-clause BSD style license - see LICENSE

import argparse
import logging
import warnings
from pathlib import Path

import romancal.orientation.set_telescope_pointing as stp

# Configure logging
logger = logging.getLogger("romancal")
logger_handler = logging.StreamHandler()
logger.addHandler(logger_handler)
logger_format_debug = logging.Formatter(
    "%(levelname)s:%(filename)s::%(funcName)s: %(message)s"
)


def main():
    """Set the initial world coordinate system."""
    parser = argparse.ArgumentParser(
        description=(
            "Update basic WCS information in Roman exposures from the engineering database."
            " For detailed information, see"
            " https://roman-pipeline.readthedocs.io/en/latest/roman/orientation/set_telescope_pointing.html"
        )
    )
    parser.add_argument(
        "exposure", type=str, nargs="+", help="List of Roman exposures to update."
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity. Specifying multiple times adds more output.",
    )
    parser.add_argument(
        "--allow-default",
        action="store_true",
        help="If pointing information cannot be determine, use header information.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform all actions but do not save the results",
    )
    parser.add_argument(
        "--save-transforms", action="store_true", help="Save transforms."
    )
    parser.add_argument(
        "--tolerance",
        type=int,
        default=60,
        help="Seconds beyond the observation time to search for telemetry. Default: %(default)s",
    )
    parser.add_argument(
        "--engdb_url",
        type=str,
        default=None,
        help=(
            "URL of the engineering database."
            " If not specified, the environment variable 'ENG_BASE_URL' is used."
            " Otherwise, a hardwired default is used."
        ),
    )

    args = parser.parse_args()

    # Set output detail.
    level = stp.LOGLEVELS[min(len(stp.LOGLEVELS) - 1, args.verbose)]
    logger.setLevel(level)
    if level <= logging.DEBUG:
        logger_handler.setFormatter(logger_format_debug)
    logger.info("set_telescope_pointing called with args %s", args)

    # Calculate WCS for all inputs.
    for filename in args.exposure:
        logger.info("")
        logger.info("------")
        logger.info(f"Setting pointing for {filename}")

        # Create path for saving the transforms.
        transform_path = None
        if args.save_transforms:
            path = Path(filename)
            transform_path = path.with_name(f"{path.stem}_transforms.asdf")

        try:
            stp.add_wcs(
                filename,
                engdb_url=args.engdb_url,
                tolerance=args.tolerance,
                allow_default=args.allow_default,
                dry_run=args.dry_run,
                save_transforms=transform_path,
            )
        except (TypeError, ValueError) as exception:
            logger.warning("Cannot determine pointing information: %s", str(exception))
            logger.debug("Full exception:", exc_info=exception)


def deprecated_name():
    """Raise warning if filename.* is no longer used, and provide correct one."""
    filename = Path(__file__)
    warnings.warn(
        f"usage of `{filename.name}` is deprecated; use `{filename.stem}` instead",
        stacklevel=2,
    )

    main()


if __name__ == "__main__":
    main()
