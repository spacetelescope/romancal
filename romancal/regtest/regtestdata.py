import json
import os
import shutil
import sys
import time
from contextlib import chdir
from pathlib import Path

import asdf

from romancal.associations._load_asn import load_asn
from romancal.regtest._differ import Diff, DiffSettings

# Define a request timeout in seconds
TIMEOUT = int(os.environ.get("TEST_BIGDATA_TIMEOUT", 30))
CHUNK_SIZE = int(os.environ.get("TEST_BIGDATA_CHUNK_SIZE", 1048576))
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
        self._bigdata_root_remote = os.environ.get("TEST_BIGDATA_REMOTE")
        self._bigdata_refresh = bool(
            json.loads(os.environ.get("TEST_BIGDATA_REFRESH", "False").lower())
        )

        if self._bigdata_refresh and not self._bigdata_root_remote:
            raise ValueError(
                "The TEST_BIGDATA_REMOTE environment variable must be defined if TEST_BIGDATA_REFRESH is enabled"
            )

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

    def _get_bigdata_remote(self, src_url, dest_url):
        _download(src_url, dest_url)

    def _get_bigdata_local(self, src_path, dest_path):
        if src_path == dest_path:
            raise OSError(f"Source and destination paths are identical: {src_path}")
        shutil.copy2(src_path, dest_path)

    def _get_bigdata(self, path):
        src = os.path.join(self._bigdata_root, self._inputs_root, self._env, path)
        filename = os.path.basename(src)
        dest = os.path.abspath(os.path.join(os.curdir, filename))

        # if src is path, if path doesn't exist OR refresh, fetch from remote

        if src.startswith("https://"):
            # when src is a url, always fetch
            self._get_bigdata_remote(src, dest)
        else:
            # conditionally fetch if refresh or input doesn't exist
            # TODO update to check if files differ?
            if self._bigdata_refresh or (
                self._bigdata_root_remote and (not os.path.exists(src))
            ):
                remote_src = os.path.join(
                    self._bigdata_root_remote, self._inputs_root, self._env, path
                )
                self._get_bigdata_remote(remote_src, src)
            self._get_bigdata_local(src, dest)
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
        af = asdf.AsdfFile(
            {
                "input": self.input,
                "output": self.output,
                "truth": self.truth,
                "input_remote": self.input_remote,
                "truth_remote": self.truth_remote,
                "remote_results_path": self.remote_results_path,
                "test_name": self.test_name,
                "traceback": self.traceback,
            }
        )
        af.write_to(path)


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
