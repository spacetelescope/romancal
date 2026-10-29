import os
import pprint
from contextlib import chdir
from pathlib import Path
from textwrap import dedent

import asdf
import roman_datamodels as rdm
from ci_watson.artifactory_helpers import (
    get_bigdata,
    get_bigdata_root,
)

from romancal.associations._load_asn import load_asn
from romancal.regtest._differ import Diff, DiffSettings

# from romancal.lib.suffix import replace_suffix

# Define location of default Artifactory API key, for Jenkins use only
ARTIFACTORY_API_KEY_FILE = "/eng/ssb2/keys/svc_rodata.key"

# Define a request timeout in seconds
TIMEOUT = 30


class RegtestData:
    """Defines data paths on Artifactory and data retrieval methods"""

    def __init__(
        self,
        env="dev",
        inputs_root="roman-pipeline",
        results_root="roman-pipeline-results/regression-tests/runs/",
        docopy=True,
        input=None,
        input_remote=None,
        output=None,
        truth=None,
        truth_remote=None,
        remote_results_path=None,
        test_name=None,
        traceback=None,
        **kwargs,
    ):
        self._env = env
        self._inputs_root = inputs_root
        self._results_root = results_root
        self._bigdata_root = get_bigdata_root()

        self.docopy = docopy

        # Initialize @property attributes
        self.input = input
        self.input_remote = input_remote
        self.output = output
        self.truth = truth
        self.truth_remote = truth_remote

        # No @properties for the following attributes
        self.remote_results_path = remote_results_path
        self.test_name = test_name
        self.traceback = traceback

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

    @property
    def bigdata_root(self):
        return self._bigdata_root

    # The methods
    def get_data(self, path=None, docopy=None):
        """Copy data from Artifactory remote resource to the CWD

        Updates self.input and self.input_remote upon completion
        """
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path
        if docopy is None:
            docopy = self.docopy
        self.input = get_bigdata(self._inputs_root, self._env, path, docopy=docopy)
        self.input_remote = os.path.join(self._inputs_root, self._env, path)

        return self.input

    def get_truth(self, path=None, docopy=None):
        """Copy truth data from Artifactory remote resource to the CWD/truth

        Updates self.truth and self.truth_remote on completion
        """
        if path is None:
            path = self.truth_remote
        else:
            self.truth_remote = path
        if docopy is None:
            docopy = self.docopy

        truth_dir = Path("truth")
        truth_dir.mkdir(exist_ok=True)
        with chdir(truth_dir):
            self.truth = get_bigdata(self._inputs_root, self._env, path, docopy=docopy)
            self.truth_remote = os.path.join(self._inputs_root, self._env, path)

        return self.truth

    def get_asn(self, path=None, docopy=True, get_members=True):
        """Copy association and association members from Artifactory remote
        resource to the CWD/truth.

        Updates self.input and self.input_remote upon completion

        Parameters
        ----------
        path: str
            The remote path

        docopy : bool
            Switch to control whether or not to copy a file
            into the test output directory when running the test.
            If you wish to open the file directly from remote
            location or just to set path to source, set this to `False`.
            Default: `True`

        get_members: bool
            If an association is the input, retrieve the members.
            Otherwise, do not.
        """
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path
        if docopy is None:
            docopy = self.docopy

        # Get the association JSON file
        self.input = get_bigdata(self._inputs_root, self._env, path, docopy=docopy)
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
                    get_bigdata(
                        self._inputs_root, self._env, fullpath, docopy=self.docopy
                    )

    def to_asdf(self, path):
        tree = eval(str(self))  # noqa: S307
        af = asdf.AsdfFile(tree=tree)
        af.write_to(path)

    @classmethod
    def open(cls, filename):
        with asdf.open(filename) as af:
            return cls(**af.tree)


class DiffResult:
    def __init__(self, diff, result, truth):
        self.diff = diff
        self.result = result
        self.truth = truth

    @property
    def identical(self):
        return not self.diff

    @staticmethod
    def _file_info(file):
        try:
            with rdm.open(file) as model:
                model_type = type(model).__name__
                crds_context = (
                    model.get("meta", {})
                    .get("ref_file", {})
                    .get("crds", {})
                    .get("context", "Unknown context")
                )
        except Exception:
            model_type = "Not a model"
            crds_context = "Unknown context"
        return model_type, crds_context

    def report(self, **kwargs):
        result_model_type, result_crds_context = self._file_info(self.result)
        truth_model_type, truth_crds_context = self._file_info(self.truth)
        title = dedent(
            f"""
            Diff report for:
                result file: {self.result}
                    model type: {result_model_type}
                    crds_context: {result_crds_context}
                truth file: {self.truth}
                    model type: {truth_model_type}
                    crds_context: {truth_crds_context}
            """
        )
        report = pprint.pformat(self.diff)

        return f"{title}\n{report}"


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

    diff_result : DiffResult
        result of the comparison
    """
    settings = DiffSettings(rtol, atol, equal_nan, ignore or [])
    with asdf.open(result) as af0:
        with asdf.config_context() as cfg:
            # Disable validation of the truth file to allow for more
            # flexible comparisons when the truth file was okified from
            # a prievious "dev" version of rad (which is not strictly
            # versioned for every change between releases).
            cfg.validate_on_read = False
            with asdf.open(truth) as af1:
                return Diff(af1.tree, af0.tree, settings)
