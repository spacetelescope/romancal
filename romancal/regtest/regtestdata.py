import os
import pprint
import shutil
import sys
import time
from contextlib import chdir
from pathlib import Path

import asdf
import astropy.table
import astropy.time
import deepdiff
import gwcs
import numpy as np
import roman_datamodels as rdm
from astropy.units import Quantity
from deepdiff.operator import BaseOperator
from gwcs.wcstools import grid_from_bounding_box

from romancal.associations._load_asn import load_asn
from romancal.regtest._differ import Diff, DiffSettings

# from romancal.lib.suffix import replace_suffix

# Define location of default Artifactory API key, for Jenkins use only
ARTIFACTORY_API_KEY_FILE = "/eng/ssb2/keys/svc_rodata.key"

# Define a request timeout in seconds
TIMEOUT = int(os.environ.get("TEST_BIGDATA_TIMEOUT", 30))
CHUNK_SIZE = int(os.environ.get("TEST_BIGDATA_CHUNK_SIZE", 16384))
RETRY_MAX = int(os.environ.get("TEST_BIGDATA_RETRY_MAX", 3))
RETRY_DELAY = int(os.environ.get("TEST_BIGDATA_RETRY_DELAY", 5))


def _download(url, dest, timeout=TIMEOUT, chunk_size=CHUNK_SIZE):
    """Simple HTTP/HTTPS downloader."""
    # Optional import: requests is not needed for local big data setup.
    import requests

    dest = os.path.abspath(dest)

    retries = RETRY_MAX

    while retries:
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with open(dest, "w+b") as data:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        data.write(chunk)
            return dest
        except Exception as err:
            retries -= 1
            if not retries:
                raise
            print(f"download: {url} failed({err}), trying again...", file=sys.stderr)
            time.sleep(RETRY_DELAY)


class RegtestData:
    """Defines data paths on Artifactory and data retrieval methods"""

    def __init__(
        self,
        env,
        inputs_root,
        results_root,
    ):
        self._env = env
        self._inputs_root = inputs_root
        self._results_root = results_root
        self._bigdata_root = os.environ.get("TEST_BIGDATA")

        # Initialize @property attributes
        self.input = None
        self.input_remote = None
        self.output = None
        self.truth = None
        self.truth_remote = None

        # No @properties for the following attributes
        self.remote_results_path = None
        self.test_name = None
        self.traceback = None

        # Initialize non-initialized attributes
        self.asn = None

    def __repr__(self):
        return pprint.pformat(
            dict(
                input=self.input,
                output=self.output,
                truth=self.truth,
                input_remote=self.input_remote,
                truth_remote=self.truth_remote,
                remote_results_path=self.remote_results_path,
                test_name=self.test_name,
                traceback=self.traceback,
            ),
            indent=1,
        )

    @property
    def input_remote(self):
        if self._input_remote is not None:
            return os.path.join(*self._input_remote)
        else:
            return None

    @input_remote.setter
    def input_remote(self, value):
        if value:
            self._input_remote = value.split(os.sep)
        else:
            self._input_remote = value

    @property
    def truth_remote(self):
        if self._truth_remote is not None:
            return os.path.join(*self._truth_remote)
        else:
            return None

    @truth_remote.setter
    def truth_remote(self, value):
        if value:
            self._truth_remote = value.split(os.sep)
        else:
            self._truth_remote = value

    @property
    def input(self):
        return self._input

    @input.setter
    def input(self, value):
        if value:
            self._input = os.path.abspath(value)
        else:
            self._input = value

    @property
    def truth(self):
        return self._truth

    @truth.setter
    def truth(self, value):
        if value:
            self._truth = os.path.abspath(value)
        else:
            self._truth = value

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, value):
        if value:
            self._output = os.path.abspath(value)
        else:
            self._output = value

    def _get_bigdata(self, path):
        src = os.path.join(self._bigdata_root, self._inputs_root, self._env, path)

        filename = os.path.basename(src)
        dest = os.path.abspath(os.path.join(os.curdir, filename))

        if os.path.exists(src):
            if src == dest:
                raise OSError(f"Source and destination paths are identical: {src}")
            shutil.copy2(src, dest)
        else:
            _download(src, dest)
        return dest

    # The methods
    def get_data(self, path=None):
        """Copy data from Artifactory remote resource to the CWD

        Updates self.input and self.input_remote upon completion
        """
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path
        self.input = self._get_bigdata(path)
        self.input_remote = os.path.join(self._inputs_root, self._env, path)

        return self.input

    def get_truth(self, path=None):
        """Copy truth data from Artifactory remote resource to the CWD/truth

        Updates self.truth and self.truth_remote on completion
        """
        if path is None:
            path = self.truth_remote
        else:
            self.truth_remote = path

        truth_dir = Path("truth")
        truth_dir.mkdir(exist_ok=True)
        with chdir(truth_dir):
            self.truth = self._get_bigdata(path)
            self.truth_remote = os.path.join(self._inputs_root, self._env, path)

        return self.truth

    def get_asn(self, path=None, get_members=True):
        """Copy association and association members from Artifactory remote
        resource to the CWD/truth.

        Updates self.input and self.input_remote upon completion

        Parameters
        ----------
        path: str
            The remote path

        get_members: bool
            If an association is the input, retrieve the members.
            Otherwise, do not.
        """
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path

        # Get the association JSON file
        self.input = self._get_bigdata(path)
        with open(self.input) as fp:
            asn = load_asn(fp)
            self.asn = asn

        # Get each member in the association as well
        if get_members:
            for product in asn["products"]:
                for member in product["members"]:
                    fullpath = os.path.join(
                        os.path.dirname(self.input_remote), member["expname"]
                    )
                    self._get_bigdata(fullpath)

    def to_asdf(self, path):
        tree = eval(str(self))  # noqa: S307
        af = asdf.AsdfFile(tree=tree)
        af.write_to(path)

    @classmethod
    def open(cls, filename):
        with asdf.open(filename) as af:
            return cls(**af.tree)


def compare_asdf(result, truth, ignore=None, rtol=1e-05, atol=1e-08, equal_nan=True):
    """
    Compare 2 asdf files: result and truth. Note that this comparison is
    asymmetric (swapping result and truth will give a different result).

    Parameters
    ----------

    result : str
        Filename of result asdf file

    truth : str
        Filename of truth asdf file

    ignore : list
        List of tree node paths to ignore during the comparison

    rtol : float
        rtol argument passed to `numpyp.allclose`

    atol : float
        atol argument passed to `numpyp.allclose` and ``DeepDiff``
        as ``math_epsilon``.

    equal_nan : bool
        Ignore nan inequality

    Returns
    -------

    diff_result : Diff
        result of the comparison
    """
    settings = DiffSettings(rtol, atol, equal_nan, set(ignore or []))
    with asdf.open(result) as af0:
        with asdf.config_context() as cfg:
            # Disable validation of the truth file to allow for more
            # flexible comparisons when the truth file was okified from
            # a prievious "dev" version of rad (which is not strictly
            # versioned for every change between releases).
            cfg.validate_on_read = False
            with asdf.open(truth) as af1:
                return Diff(af1.tree, af0.tree, settings)
