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
        pattern = re.compile(
            r".*_(?P<skycells>[0-9p]*x[0-9]*y[0-9]*)_f[0-9]*_coadd\.asdf$"
        )
        groups = {}
        for filename in filelist:
            match = pattern.match(filename)
            if match:
                key = match.group("skycells")
                groups.setdefault(key, []).append(filename)
        return groups

    def create_multiband_asn(self):
        """
        Create a multiband association from a list of files.

        Parameters:
        files (list): List of file paths or pattern to include in the association.

        Returns:
        dict: The created association.
        """
        for skycell_id, filenames in self.skycell_groups.items():
            # Get prefixes for all combinations of data_release_id and product_type from filenames
            # (r00001_{data_release_id}_{product_type}_{skycell_id}_asn.json)
            prefixes = {x.split(skycell_id)[0] for x in filenames}
            for prefix in prefixes:
                # Get all files that match this prefix (data_release_id + product_type) and skycell
                files = [x for x in filenames if x.startswith(f"{prefix}{skycell_id}")]
                args = [
                    *files,
                    "-o",
                    f"{prefix}{skycell_id}_asn.json",
                    "--product-name",
                    f"{skycell_id}",
                    "--data-release-id",
                    prefix.split("_")[1],
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
