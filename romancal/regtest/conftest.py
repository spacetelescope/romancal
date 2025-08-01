import copy
import getpass
import json
import os
from datetime import datetime

import pytest
from ci_watson.artifactory_helpers import UPLOAD_SCHEMA

from romancal.regtest.regtestdata import RegtestData
from romancal.regtest.resource_tracker import ResourceTracker

TODAYS_DATE = datetime.now().strftime("%Y-%m-%d")


@pytest.fixture(scope="session")
def artifactory_repos(pytestconfig):
    """Provides Artifactory inputs_root and results_root"""
    inputs_root = pytestconfig.getini("inputs_root")
    if not inputs_root:
        inputs_root = "roman-pipeline"
    else:
        inputs_root = inputs_root[0]

    results_root = pytestconfig.getini("results_root")
    if not results_root:
        results_root = "roman-pipeline-results/regression-tests/runs/"
    else:
        results_root = results_root[0]
    return inputs_root, results_root


@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Create request.node.report_setup and request.node.report_call to be
    used in generate_artifactory_json() fixture.
    """
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()

    # set a report attribute for each phase of a call, which can
    # be "setup", "call", "teardown"
    setattr(item, "report_" + rep.when, rep)


def postmortem(request, fixturename):
    """Retrieve a fixture object if a test failed, else return None"""
    try:
        if request.node.report_setup.passed:
            try:
                if request.node.report_call.failed:
                    return request.node.funcargs.get(fixturename, None)
            # Handle case where `report_call` hasn't been generated yet
            # because the test hasn't finished running, as in the case of
            # a user interrupt
            except AttributeError:
                return None
    except AttributeError:
        return None


def pytest_collection_modifyitems(config, items):
    # add the index of each item in the list of items
    # this is used below for artifactory_result_path
    # to produce a unique result subdirectory for
    # each test (even if that test shares a name with
    # another test which is the case for parametrized tests).
    for i, item in enumerate(items):
        item.index = i


@pytest.fixture(scope="function", autouse=True)
def generate_artifactory_json(request, artifactory_repos):
    """Pytest fixture that leaves behind JSON upload and okify specfiles
    if the rtdata or rtdata_module fixtures are used in either a test or a
    module-scoped fixture that runs a pipeline and provides the results to a
    series of test.
    """
    inputs_root, results_root = artifactory_repos

    def artifactory_result_path():
        # Generate the Artifactory result path
        whoami = getpass.getuser() or "nobody"
        user_tag = f"NOT_CI_{whoami}"
        build_tag = os.environ.get("BUILD_TAG", user_tag)
        build_matrix_suffix = os.environ.get("BUILD_MATRIX_SUFFIX", "0")
        subdir = f"{TODAYS_DATE}_{build_tag}_{build_matrix_suffix}"
        testname = request.node.originalname or request.node.name
        basename = f"{request.node.index}_{testname}"

        return os.path.join(results_root, subdir, basename) + os.sep

    yield
    # Execute the following at test teardown
    upload_schema_pattern = []
    okify_schema_pattern = []

    rtdata = postmortem(request, "rtdata") or postmortem(request, "rtdata_module")
    if rtdata:
        try:
            # The function_jail fixture uses tmp_path
            cwd = str(request.node.funcargs["tmp_path"])
        except KeyError:
            # The module_jail fixture (module-scoped) returns the path
            cwd = str(request.node.funcargs["module_jail"])
        rtdata.remote_results_path = artifactory_result_path()
        rtdata.test_name = request.node.name
        # Dump the failed test traceback into rtdata
        rtdata.traceback = str(request.node.report_call.longrepr)

        # Upload and allow okify of truth by rtdata.output, if the test did not
        # fail before producing rtdata.output
        if rtdata.output and os.path.exists(rtdata.output):
            # Write the rtdata class out as an ASDF file
            path_asdf = os.path.join(cwd, f"{request.node.name}_rtdata.asdf")
            rtdata.to_asdf(path_asdf)

            # Generate an OKify JSON file
            pattern = os.path.join(
                rtdata.remote_results_path, os.path.basename(rtdata.output)
            )
            okify_schema_pattern.append(pattern)
            okify_schema = generate_upload_schema(
                okify_schema_pattern, f"{os.path.dirname(rtdata.truth_remote)}/"
            )

            jsonfile = os.path.join(cwd, f"{request.node.name}_okify.json")
            with open(jsonfile, "w") as fd:
                json.dump(okify_schema, fd, indent=2)

            # Generate an upload JSON file, including the OKify, asdf file
            upload_schema_pattern.append(rtdata.output)
            upload_schema_pattern.append(os.path.abspath(jsonfile))
            upload_schema_pattern.append(path_asdf)
            upload_schema = generate_upload_schema(
                upload_schema_pattern, rtdata.remote_results_path
            )

            jsonfile = os.path.join(cwd, f"{request.node.name}_results.json")
            with open(jsonfile, "w") as fd:
                json.dump(upload_schema, fd, indent=2)


def generate_upload_schema(pattern, target, recursive=False):
    """
    Generate JSON schema for Artifactory upload specfile using JFROG.

    Docs can be found at
    https://www.jfrog.com/confluence/display/RTF/Using+File+Specs

    Parameters
    ----------
    pattern : str or list of strings
        Specifies the local file system path to test results which should be
        uploaded to Artifactory. You can specify multiple artifacts by using
        wildcards or a regular expression as designated by the regexp property.

    target : str
        Specifies the target path in Artifactory in the following format::

            [repository_name]/[repository_path]

    recursive : bool, optional
        Specify whether or not to identify files listed in sub-directories
        for uploading.  Default: `False`

    Returns
    -------
    upload_schema : dict
        Dictionary specifying the upload schema
    """
    recursive = repr(recursive).lower()

    if not isinstance(pattern, str):
        # Populate schema for this test's data
        upload_schema = {"files": []}

        for p in pattern:
            temp_schema = copy.deepcopy(UPLOAD_SCHEMA["files"][0])
            temp_schema.update({"pattern": p, "target": target, "recursive": recursive})
            upload_schema["files"].append(temp_schema)
    else:
        # Populate schema for this test's data
        upload_schema = copy.deepcopy(UPLOAD_SCHEMA)
        upload_schema["files"][0].update(
            {"pattern": pattern, "target": target, "recursive": recursive}
        )
    return upload_schema


def _rtdata_fixture_implementation(artifactory_repos, envopt, request):
    """Provides the RemoteResource class"""
    inputs_root, results_root = artifactory_repos
    rtdata = RegtestData(env=envopt, inputs_root=inputs_root, results_root=results_root)

    yield rtdata


@pytest.fixture(scope="function")
def rtdata(artifactory_repos, envopt, request, function_jail):
    yield from _rtdata_fixture_implementation(artifactory_repos, envopt, request)


@pytest.fixture(scope="module")
def rtdata_module(artifactory_repos, envopt, request, module_jail):
    yield from _rtdata_fixture_implementation(artifactory_repos, envopt, request)


@pytest.fixture
def ignore_asdf_paths():
    ignore_attr = [
        # generic asdf stuff that will contain program version numbers
        # and other things that will almost certainly change in every case
        "asdf_library",
        "history",
        # roman-specific stuff to ignore
        "roman.meta.ref_file.crds.version",
        "roman.meta.calibration_software_version",
        "roman.cal_logs",
        "roman.meta.cal_logs",
        "roman.meta.date",
        "roman.meta.file_date",
        "roman.individual_image_cal_logs",
        "roman.meta.individual_image_meta",
    ]

    return {"ignore": ignore_attr}


@pytest.fixture(scope="module")
def resource_tracker():
    """Fixture to return the current module-scoped ResourceTracker.

    Use by calling ``track`` to generate a context in which resource
    usage will be tracked.

    >>>
    >> with resource_tracker.track():
    >>     # do stuff

    For resources used during tests providing a function-scoped
    request fixture result as the log argument will also log the
    used resources to the junit results.xml.

    >>>
    >> def test_something(resource_tracker, request):
    >>     with resource_tracker.track(log=request):
    >>         # do stuff

    For resources used during fixtures the tracked resources
    can be logged in a separate test using ``log_tracked_resources``.
    """
    return ResourceTracker()


@pytest.fixture()
def log_tracked_resources(resource_tracker, request):
    """Fixture to log resources tracked by ``resource_tracker``.

    >>>
    >> @pytest.fixture
    >> def my_fixture(resource_tracker):
    >>     with resource_tracker.track():
    >>         # do stuff
    >>
    >> def test_write_log(log_tracked_resources, my_fixture):
    >>     log_tracked_resources()
    """

    def callback():
        resource_tracker.log(request)

    yield callback
