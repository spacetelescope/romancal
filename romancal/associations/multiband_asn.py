import argparse
import glob
import logging
import re

from . import asn_from_list

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel("INFO")

__all__ = ["MultibandAssociation"]


_COADD_RE = re.compile(
    r"^(?P<prefix>.*_)(?P<skycell>[0-9p]*x[0-9]*y[0-9]*)_(?P<filter>f[0-9]+)_coadd\.asdf$"
)


class MultibandAssociation:
    """A class to create multiband associations."""

    def __init__(self, files):
        self.files = self._parse_file_list(files)
        self.skycell_groups = self._get_skycell_groups(self.files)

    def _parse_file_list(self, files):
        """
        Parse a file list, expanding wildcards if present.

        If the input list contains a single string with wildcard characters
        ('*' or '?'), expand it to the matching files using glob. Otherwise,
        return the list as is.

        Parameters
        ----------
        files : list of str
            List of file paths or a single wildcard pattern.

        Returns
        -------
        list of str
            List of file paths, expanded if a wildcard was provided.
        """
        if len(files) == 1 and any(char in files[0] for char in ["*", "?"]):
            return glob.glob(files[0])
        return files

    def _get_skycell_groups(self, filelist):
        """
        Create skycell groups based on the unique skycell identifiers from a list of filenames.

        Parameters
        ----------
        filelist : list of str
            List of filenames.

        Returns
        -------
        dict
            Dictionary mapping skycell identifiers to lists of filenames.
        """
        groups = {}
        for filename in filelist:
            match = _COADD_RE.match(filename)
            if match:
                key = match.group("skycell")
                groups.setdefault(key, []).append(filename)
        return groups

    def create_multiband_asn(self):
        """Create a multiband association from a list of files."""
        for skycell_id, filenames in self.skycell_groups.items():
            # Group by the file prefix (everything up to and including the underscore
            # before the skycell id). This prefix encodes program, data_release_id,
            # and product type (e.g. visit/full/pass).
            prefixes = set()
            for filename in filenames:
                match = _COADD_RE.match(filename)
                if match:
                    prefixes.add(match.group("prefix"))

            for prefix in prefixes:
                # Get all files that match this prefix (data_release_id + product_type)
                # and this skycell.
                files = [x for x in filenames if x.startswith(f"{prefix}{skycell_id}")]

                # Determine data release id from the prefix: rPPPPP_<data_release_id>_...
                # If parsing fails, fall back to 'p'.
                try:
                    data_release_id = prefix.split("_")[1]
                except IndexError:
                    data_release_id = "p"

                # Roman naming conventions for the archive catalog (Level 4):
                # - Multiband catalogs always omit the optical element (filter)
                #   regardless of whether they contain single or multiple filters.
                # - This distinguishes multiband products from single-band products.
                product_name = f"{prefix}{skycell_id}"
                output_asn = f"{product_name}_asn.json"

                args = [
                    *files,
                    "-o",
                    output_asn,
                    "--product-name",
                    product_name,
                    "--data-release-id",
                    data_release_id,
                ]
                asn_from_list._cli(args)


def _cli():
    parser = argparse.ArgumentParser(
        description="Create a multiband association from a list of files",
        usage="multiband_asn file1.asdf file2.asdf ... fileN.asdf",
    )
    parser.add_argument(
        "files",
        type=str,
        nargs="+",
        help="List of files to include in the multiband association",
    )

    args = parser.parse_args()

    multiband_asn = MultibandAssociation(args.files)

    multiband_asn.create_multiband_asn()

    logger.info("Multiband association creation complete.")
