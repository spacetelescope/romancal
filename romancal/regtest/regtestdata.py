import os
import os.path as op
import pprint
import shutil
import sys
from difflib import unified_diff
from glob import glob as _sys_glob
from pathlib import Path
from textwrap import dedent

import asdf
import astropy.table
import astropy.time
import deepdiff
import gwcs
import numpy as np
import requests
from astropy.units import Quantity
from ci_watson.artifactory_helpers import (
    BigdataError,
    check_url,
    get_bigdata,
    get_bigdata_root,
)
from deepdiff.operator import BaseOperator
from gwcs.wcstools import grid_from_bounding_box

# from romancal.lib.suffix import replace_suffix
from romancal.stpipe import RomanStep

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
        results_root="roman-pipeline-results",
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

    def data_glob(self, path=None, glob="*", docopy=None):
        """Get a list of files"""
        if path is None:
            path = self.input_remote
        else:
            self.input_remote = path
        if docopy is None:
            docopy = self.docopy

        # Get full path and proceed depending on whether
        # is a local path or URL.
        root = self.bigdata_root
        if op.exists(root):
            root_path = op.join(root, self._inputs_root, self._env)
            root_len = len(root_path) + 1
            path = op.join(root_path, path)
            file_paths = _data_glob_local(path, glob)
        elif check_url(root):
            root_len = len(self._env) + 1
            file_paths = _data_glob_url(
                self._inputs_root, self._env, path, glob, root=root
            )
        else:
            raise BigdataError(f"Path cannot be found: {path}")

        # Remove the root from the paths
        file_paths = [file_path[root_len:] for file_path in file_paths]
        return file_paths

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
        os.makedirs("truth", exist_ok=True)
        os.chdir("truth")
        try:
            self.truth = get_bigdata(self._inputs_root, self._env, path, docopy=docopy)
            self.truth_remote = os.path.join(self._inputs_root, self._env, path)
        except BigdataError:
            os.chdir("..")
            raise
        os.chdir("..")

        return self.truth

    # def get_asn(self, path=None, docopy=True, get_members=True):
    #     """Copy association and association members from Artifactory remote
    #     resource to the CWD/truth.
    #
    #     Updates self.input and self.input_remote upon completion
    #
    #     Parameters
    #     ----------
    #     path: str
    #         The remote path
    #
    #     docopy : bool
    #         Switch to control whether or not to copy a file
    #         into the test output directory when running the test.
    #         If you wish to open the file directly from remote
    #         location or just to set path to source, set this to `False`.
    #         Default: `True`
    #
    #     get_members: bool
    #         If an association is the input, retrieve the members.
    #         Otherwise, do not.
    #     """
    #     if path is None:
    #         path = self.input_remote
    #     else:
    #         self.input_remote = path
    #     if docopy is None:
    #         docopy = self.docopy
    #
    #     # Get the association JSON file
    #     self.input = get_bigdata(self._inputs_root, self._env, path,
    #         docopy=docopy)
    #     with open(self.input) as fp:
    #         asn = load_asn(fp)
    #         self.asn = asn
    #
    #     # Get each member in the association as well
    #     if get_members:
    #         for product in asn['products']:
    #             for member in product['members']:
    #                 fullpath = os.path.join(
    #                     os.path.dirname(self.input_remote),
    #                     member['expname'])
    #                 get_bigdata(self._inputs_root, self._env, fullpath,
    #                             docopy=self.docopy)

    def to_asdf(self, path):
        tree = eval(str(self))
        af = asdf.AsdfFile(tree=tree)
        af.write_to(path)

    @classmethod
    def open(cls, filename):
        with asdf.open(filename) as af:
            return cls(**af.tree)


def run_step_from_dict(rtdata, **step_params):
    """Run Steps with given parameter

    Parameters
    ----------
    rtdata: RegtestData
        The artifactory instance

    step_params: dict
        The parameters defining what step to run with what input

    Returns
    -------
    rtdata: RegtestData
        Updated `RegtestData` object with inputs set.

    Notes
    -----
    `step_params` looks like this:
    {
        'input_path': str or None  # The input file path, relative to artifactory
        'step': str                # The step to run, either a class or a config file
        'args': list,              # The arguments passed to `Step.from_cmdline`
    }
    """

    # Get the data. If `step_params['input_path]` is not
    # specified, the presumption is that `rtdata.input` has
    # already been retrieved.
    # input_path = step_params.get('input_path', None)
    # if input_path:
    #     try:
    #         rtdata.get_asn(input_path)
    #     except AssociationNotValidError:
    #         rtdata.get_data(input_path)

    # Figure out whether we have a config or class
    step = step_params["step"]
    if step.endswith((".asdf", ".cfg")):
        step = os.path.join("config", step)

    # Run the step
    full_args = [step, rtdata.input]
    full_args.extend(step_params["args"])

    RomanStep.from_cmdline(full_args)

    return rtdata


def run_step_from_dict_mock(rtdata, source, **step_params):
    """Pretend to run Steps with given parameter but just copy data

    For long running steps where the result already exists, just
    copy the data from source

    Parameters
    ----------
    rtdata: RegtestData
        The artifactory instance

    step_params: dict
        The parameters defining what step to run with what input

    source: Path-like folder
        The folder to copy from. All regular files are copied.

    Returns
    -------
    rtdata: RegtestData
        Updated `RegtestData` object with inputs set.

    Notes
    -----
    `step_params` looks like this:
    {
        'input_path': str or None  # The input file path, relative to artifactory
        'step': str                # The step to run, either a class or a config file
        'args': list,              # The arguments passed to `Step.from_cmdline`
    }
    """

    # Get the data. If `step_params['input_path]` is not
    # specified, the presumption is that `rtdata.input` has
    # already been retrieved.
    # input_path = step_params.get('input_path', None)
    # if input_path:
    #     try:
    #         rtdata.get_asn(input_path)
    #     except AssociationNotValidError:
    #         rtdata.get_data(input_path)

    # Copy the data
    for file_name in os.listdir(source):
        file_path = os.path.join(source, file_name)
        if os.path.isfile(file_path):
            shutil.copy(file_path, ".")

    return rtdata


def is_like_truth(rtdata, ignore_asdf_paths, output, truth_path, is_suffix=True):
    """Compare step outputs with truth

    Parameters
    ----------
    rtdata: RegtestData
        The artifactory object from the step run.

    ignore_asdf_paths: dict
        The asdf `diff` keyword arguments

    output: str
        The suffix or full file name to check on.

    truth_path: str
        Location of the truth files.

    is_suffix: bool
        Interpret `output` as just a suffix on the expected output root.
        Otherwise, assume it is a full file name
    """
    __tracebackhide__ = True
    # If given only a suffix, get the root to change the suffix of.
    # If the input was an association, the output should be the name of
    # the product. Otherwise, output is based on input.
    if is_suffix:
        # suffix = output
        if rtdata.asn:
            output = rtdata.asn["products"][0]["name"]
        else:
            output = os.path.splitext(os.path.basename(rtdata.input))[0]
        # output = replace_suffix(output, suffix) + '.asdf'
    rtdata.output = output

    rtdata.get_truth(os.path.join(truth_path, output))

    # diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    report = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert report is None, report


def text_diff(from_path, to_path):
    """Assertion helper for diffing two text files

    Parameters
    ----------
    from_path: str
        File to diff from.

    to_path: str
        File to diff to.  The truth.

    Returns
    -------
    diffs: [str[,...]]
        A generator of a list of strings that are the differences.
        The output from `difflib.unified_diff`
    """
    __tracebackhide__ = True
    with open(from_path) as fh:
        from_lines = fh.readlines()
    with open(to_path) as fh:
        to_lines = fh.readlines()

    diffs = unified_diff(from_lines, to_lines, from_path, to_path)

    diff = list(diffs)
    if len(diff) > 0:
        diff.insert(0, "\n")
        raise AssertionError("".join(diff))
    else:
        return True


def _data_glob_local(*glob_parts):
    """Perform a glob on the local path

    Parameters
    ----------
    glob_parts: (path-like,[...])
        List of components that will be built into a single path

    Returns
    -------
    file_paths: [str[, ...]]
        Full file paths that match the glob criterion
    """
    full_glob = Path().joinpath(*glob_parts)
    return _sys_glob(str(full_glob))


def _data_glob_url(*url_parts, root=None):
    """
    Parameters
    ----------
    url: (str[,...])
        List of components that will be used to create a URL path

    root: str
        The root server path to the Artifactory server.
        Normally retrieved from `get_bigdata_root`.

    Returns
    -------
    url_paths: [str[, ...]]
        Full URLS that match the glob criterion
    """
    # Fix root root-ed-ness
    if root.endswith("/"):
        root = root[:-1]

    # Access
    try:
        envkey = os.environ["API_KEY_FILE"]
    except KeyError:
        envkey = ARTIFACTORY_API_KEY_FILE

    try:
        with open(envkey) as fp:
            headers = {"X-JFrog-Art-Api": fp.readline().strip()}
    except (PermissionError, FileNotFoundError):
        print(
            "Warning: Anonymous Artifactory search requests are limited to "
            "1000 results. Use an API key and define API_KEY_FILE environment"
            "variable to get full search results.",
            file=sys.stderr,
        )
        headers = None

    search_url = "/".join([root, "api/search/pattern"])

    # Join and re-split the url so that every component is identified.
    url = "/".join([root] + [idx for idx in url_parts])
    all_parts = url.split("/")

    # Pick out "roman-pipeline", the repo name
    repo = all_parts[4]

    # Format the pattern
    pattern = repo + ":" + "/".join(all_parts[5:])

    # Make the query
    params = {"pattern": pattern}
    with requests.get(search_url, params=params, headers=headers, timeout=TIMEOUT) as r:
        url_paths = r.json()["files"]

    return url_paths


class NDArrayTypeOperator(BaseOperator):
    def __init__(self, rtol=1e-05, atol=1e-08, equal_nan=True, **kwargs):
        super().__init__(**kwargs)
        self.rtol = rtol
        self.atol = atol
        self.equal_nan = equal_nan

    def _compare_arrays(self, a, b):
        difference = {}
        if a.shape != b.shape:
            difference["shapes"] = [a.shape, b.shape]
        if a.dtype != b.dtype:
            difference["dtypes"] = [a.dtype, b.dtype]
        if isinstance(a, Quantity) and isinstance(b, Quantity):
            if a.unit != b.unit:
                difference["units"] = [a.unit, b.unit]
        if a.dtype.names or b.dtype.names:
            if a.dtype.names != b.dtype.names:
                difference["names"] = [a.dtype.names, b.dtype.names]
            if not difference:
                # only compare values if no other difference exists
                column_differences = {}
                for name in a.dtype.names:
                    column_difference = self._compare_arrays(a[name], b[name])
                    if column_difference:
                        column_differences[name] = column_difference
                if column_differences:
                    difference["column_values"] = column_differences
            return difference
        if not difference:  # only compare if shapes and dtypes match
            if not np.allclose(
                a, b, rtol=self.rtol, atol=self.atol, equal_nan=self.equal_nan
            ):
                abs_diff = np.abs(a - b)
                index = np.unravel_index(np.nanargmax(abs_diff), a.shape)
                difference["worst_abs_diff"] = {
                    "index": index,
                    "value": abs_diff[index],
                }
                with np.errstate(invalid="ignore", divide="ignore"):
                    # 0 / 0 == nan and produces an 'invalid' error
                    # 1 / 0 == inf and produces a 'divide' error
                    # ignore these here for computing the fractional diff
                    fractional_diff = np.abs(a / b)
                index = np.unravel_index(np.nanargmax(fractional_diff), a.shape)
                difference["worst_fractional_diff"] = {
                    "index": index,
                    "value": fractional_diff[index],
                }
                difference["abs_diff"] = np.nansum(np.abs(a - b))
                difference["n_diffs"] = np.count_nonzero(
                    np.isclose(
                        a,
                        b,
                        rtol=self.rtol,
                        atol=self.atol,
                        equal_nan=self.equal_nan,
                    )
                )
        return difference

    def give_up_diffing(self, level, diff_instance):
        difference = self._compare_arrays(level.t1, level.t2)
        if difference:
            diff_instance.custom_report_result("arrays_differ", level, difference)
        return True


class TableOperator(NDArrayTypeOperator):
    def give_up_diffing(self, level, diff_instance):
        a = level.t1.as_array()
        b = level.t2.as_array()
        difference = self._compare_arrays(a, b)
        # also compare meta  for tables
        if level.t1.meta != level.t2.meta:
            meta_difference = deepdiff.DeepDiff(
                level.t1.meta,
                level.t2.meta,
                ignore_nan_inequality=self.equal_nan,
                math_epsilon=self.atol,
            )
            if meta_difference:
                difference["metas_differ"] = meta_difference
        if difference:
            diff_instance.custom_report_result("tables_differ", level, difference)
        return True


class TimeOperator(BaseOperator):
    def give_up_diffing(self, level, diff_instance):
        if level.t1 != level.t2:
            diff_instance.custom_report_result(
                "times_differ",
                level,
                {
                    "difference": level.t1 - level.t2,
                },
            )
        return True


def _wcs_to_ra_dec(wcs):
    x, y = grid_from_bounding_box(wcs.bounding_box)
    return wcs(x, y)


class WCSOperator(NDArrayTypeOperator):
    def give_up_diffing(self, level, diff_instance):
        # for comparing wcs instances this function evaluates
        # each wcs and compares the resulting ra and dec outputs
        # TODO should we compare the bounding boxes?
        ra_a, dec_a = _wcs_to_ra_dec(level.t1)
        ra_b, dec_b = _wcs_to_ra_dec(level.t2)
        ra_diff = self._compare_arrays(ra_a, ra_b)
        dec_diff = self._compare_arrays(dec_a, dec_b)
        difference = {}
        if ra_diff:
            difference["ra"] = ra_diff
        if dec_diff:
            difference["dec"] = dec_diff
        if difference:
            diff_instance.custom_report_result("wcs_differ", level, difference)
        return True


class DiffResult:
    def __init__(self, diff, result, truth):
        self.diff = diff
        self.result = result
        self.truth = truth

    @property
    def identical(self):
        return not self.diff

    @staticmethod
    def _model_type(file):
        import roman_datamodels as rdm

        try:
            with rdm.open(file) as mdl:
                return mdl.__class__.__name__
        except Exception:
            return "Not a model"

    def report(self, **kwargs):
        title = dedent(
            f"""
            Diff report for:
                result file: {self.result}
                    model type: {self._model_type(self.result)}
                truth file: {self.truth}
                    model type: {self._model_type(self.result)}
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
    exclude_paths = []
    ignore = ignore or []
    for path in ignore:
        key_path = "".join([f"['{k}']" for k in path.split(".")])
        exclude_paths.append(f"root{key_path}")
    operators = [
        NDArrayTypeOperator(
            rtol,
            atol,
            equal_nan,
            types=[asdf.tags.core.NDArrayType, np.ndarray],
        ),
        TimeOperator(types=[astropy.time.Time]),
        TableOperator(rtol, atol, equal_nan, types=[astropy.table.Table]),
        WCSOperator(rtol, atol, equal_nan, types=[gwcs.WCS]),
    ]
    # warnings can be seen in regtest runs which indicate
    # that ddtrace logs are evaluated at times after the below
    # with statement exits resulting in access attempts on the
    # closed asdf file. To try and avoid that we disable
    # lazy loading and memmory mapping
    open_kwargs = {
        "lazy_load": False,
        "copy_arrays": True,
    }
    with (
        asdf.open(result, **open_kwargs) as af0,
        asdf.open(truth, **open_kwargs) as af1,
    ):
        # swap the inputs here so DeepDiff(truth, result)
        # this will create output with 'new_value' referring to
        # the value in the result and 'old_value' referring to the truth
        diff = deepdiff.DeepDiff(
            af1.tree,
            af0.tree,
            ignore_nan_inequality=equal_nan,
            custom_operators=operators,
            exclude_paths=exclude_paths,
            math_epsilon=atol,
            ignore_type_in_groups=[asdf.tags.core.NDArrayType, np.ndarray],
        )
        return DiffResult(diff, result, truth)
