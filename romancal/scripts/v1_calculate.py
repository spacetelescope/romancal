#!/usr/bin/env python

"""V1 Calculation based on time and engineering database info."""

import argparse
import logging
from sys import stdout

from astropy.time import Time

import romancal.orientation.set_telescope_pointing as stp
from romancal.lib.engdb.engdb_tools import AVAILABLE_SERVICES
from romancal.orientation import v1_calculate

# Available reduce functions
REDUCE_FUNCS_MAPPING = {
    "all": stp.all_pointings,
    "first": stp.first_pointing,
    "average": stp.pointing_from_average,
}
REDUCE_FUNCS = list(REDUCE_FUNCS_MAPPING.keys())


# Begin execution
def main():
    """Calculate V1, the telescope's boresight axis, the direction the telescope is pointing."""
    parser = argparse.ArgumentParser(description="Calculate V1 over a time period.")

    parser.add_argument(
        "time_sources",
        type=str,
        nargs="+",
        help=(
            "Either a list of Roman data files to retrieve the timing from"
            " or a start and end time to retrieve pointing information for."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=stdout,
        help="File to write V1 calculation table to. Default is standard output.",
    )
    parser.add_argument(
        "--pointing",
        type=str,
        choices=REDUCE_FUNCS,
        default="average",
        help=(
            "Which pointing(s) to use within the specified time interval."
            f" Choices: {REDUCE_FUNCS}"
            " (default: %(default)s)"
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity. Specifying multiple times adds more output.",
    )
    parser.add_argument(
        "--service",
        type=str,
        default="mast",
        choices=[name for name in AVAILABLE_SERVICES],
        help="Database service to use. Default: %(default)s",
    )

    # Arguments pertinent only to the EngdbMast service.
    if "mast" in AVAILABLE_SERVICES:
        parser.add_argument(
            "--engdb-url",
            type=str,
            default=None,
            help=(
                "URL of the engineering database."
                " If not specified, the environment variable 'ENG_BASE_URL' is used."
                " Otherwise, a hardwired default is used."
            ),
        )

    # Arguments pertinent only to the EngdbEDP service
    if "edp" in AVAILABLE_SERVICES:
        parser.add_argument(
            "--environment",
            type=str,
            default="test",
            choices=["dev", "test", "int", "ops"],
            help="Operational environment in use. Default: %(default)s",
        )
        parser.add_argument(
            "--path-to-cc",
            type=str,
            help="Full path to the required kerberos authentication keytab file",
        )

    args = parser.parse_args()

    # Set output detail.
    logger = logging.getLogger("romancal")
    logger_handler = logging.StreamHandler()
    logger.addHandler(logger_handler)
    logger_format_debug = logging.Formatter(
        "%(levelname)s:%(filename)s::%(funcName)s: %(message)s"
    )
    level = stp.LOGLEVELS[min(len(stp.LOGLEVELS) - 1, args.verbose)]
    logger.setLevel(level)
    if level <= logging.DEBUG:
        logger_handler.setFormatter(logger_format_debug)

    # Gather the service-specific args
    service_kwargs = {"service": args.service}
    for arg in ["engdb_url", "environment", "path_to_cc"]:
        try:
            service_kwargs[arg] = getattr(args, arg)
        except AttributeError:
            pass

    # Determine whether the sources are time specifications or a file list.
    if len(args.time_sources) == 2:
        try:
            obsstart = Time(args.time_sources[0])
            obsend = Time(args.time_sources[1])
        except ValueError:
            input_as_files = True
        else:
            logger.info(
                f"Retrieving V1 over the time span {obsstart.isot} - {obsend.isot}"
            )
            input_as_files = False
            if args.pointing != "all":
                logger.warning(
                    "V1 pointings have been requested over a time range. "
                    "However, the 'pointing' option is not 'all'."
                )
                logger.warning(
                    "There will only be a single result returned. Is this what was desired?"
                )
                logger.warning("Suggestion: Use '--pointing=all'")
    else:
        input_as_files = True

    # Process the file list.
    logger.info("Starting V1 calculation...")
    if input_as_files:
        v1s = v1_calculate.v1_calculate_from_models(
            args.time_sources,
            reduce_func=REDUCE_FUNCS_MAPPING[args.pointing],
            service_kwargs=service_kwargs,
        )
    else:
        v1s = v1_calculate.v1_calculate_over_time(
            obsstart,
            obsend,
            reduce_func=REDUCE_FUNCS_MAPPING[args.pointing],
            service_kwargs=service_kwargs,
        )

    formatted = v1_calculate.simplify_table(v1s)
    formatted.write(args.output, format="ascii.ecsv")
    logger.info("...V1 calculation completed.")


if __name__ == "__main__":
    main()
