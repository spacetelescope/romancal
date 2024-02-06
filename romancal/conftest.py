"""Project default for pytest"""

import inspect
import json
import os
import tempfile
from io import StringIO

import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import Shift
from gwcs import coordinate_frames as cf
from gwcs import wcs
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.assign_wcs import pointing


@pytest.fixture
def mk_tmp_dirs():
    """Create a set of temporary directorys and change to one of them."""
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()
    try:
        os.chdir(tmp_current_path)
        yield (tmp_current_path, tmp_data_path, tmp_config_path)
    finally:
        os.chdir(old_path)


@pytest.fixture
def slow(request):
    """Setup slow fixture for tests to identify if --slow
    has been specified
    """
    return request.config.getoption("--slow")


@pytest.fixture(scope="module")
def jail(request, tmpdir_factory):
    """Run test in a pristine temporary working directory, scoped to module.

    This fixture is the same as _jail in ci_watson, but scoped to module
    instead of function.  This allows a fixture using it to produce files in a
    temporary directory, and then have the tests access them.
    """
    old_dir = os.getcwd()
    path = request.module.__name__.split(".")[-1]
    if request._parent_request.fixturename is not None:
        path = path + "_" + request._parent_request.fixturename
    newpath = tmpdir_factory.mktemp(path)
    os.chdir(str(newpath))
    yield newpath
    os.chdir(old_dir)


@pytest.hookimpl(trylast=True)
def pytest_configure(config):
    terminal_reporter = config.pluginmanager.getplugin("terminalreporter")
    config.pluginmanager.register(
        TestDescriptionPlugin(terminal_reporter), "testdescription"
    )


class TestDescriptionPlugin:
    """Pytest plugin to print the test docstring when `pytest -vv` is used.

    This plug-in was added to support JWST instrument team testing and
    reporting for the JWST calibration pipeline.
    """

    def __init__(self, terminal_reporter):
        self.terminal_reporter = terminal_reporter
        self.desc = None

    def pytest_runtest_protocol(self, item):
        try:
            # Get the docstring for the test
            self.desc = inspect.getdoc(item.obj)
        except AttributeError:
            self.desc = None

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_logstart(self, nodeid, location):
        # When run as `pytest` or `pytest -v`, no change in behavior
        if self.terminal_reporter.verbosity <= 1:
            yield
        # When run as `pytest -vv`, `pytest -vvv`, etc, print the test docstring
        else:
            self.terminal_reporter.write("\n")
            yield
            if self.desc:
                self.terminal_reporter.write(f"\n{self.desc} ")


@pytest.fixture(scope="function")
def create_mock_asn_file():
    def _create_asn_file(tmp_path, members_mapping=None):
        asn_content = """
            {
                "asn_type": "None",
                "asn_rule": "DMS_ELPP_Base",
                "version_id": null,
                "code_version": "0.9.1.dev28+ge987cc9.d20230106",
                "degraded_status": "No known degraded exposures in association.",
                "program": "noprogram",
                "constraints": "No constraints",
                "asn_id": "a3001",
                "target": "none",
                "asn_pool": "test_pool_name",
                "products": [
                    {
                        "name": "files.asdf",
                        "members": [
                            {
                                "expname": "img_1.asdf",
                                "exptype": "science"
                            },
                            {
                                "expname": "img_2.asdf",
                                "exptype": "science"
                            }
                        ]
                    }
                ]
            }
        """
        if members_mapping is not None:
            asn_dict = json.loads(asn_content)
            asn_dict["products"][0]["members"] = []
            for x in members_mapping:
                asn_dict["products"][0]["members"].append(x)
            asn_content = json.dumps(asn_dict)

        asn_file_path = str(tmp_path / "sample_asn.json")
        asn_file = StringIO()
        asn_file.write(asn_content)
        with open(asn_file_path, mode="w") as f:
            print(asn_file.getvalue(), file=f)

        return asn_file_path

    return _create_asn_file


def _create_wcs(input_dm, shift_1=0, shift_2=0):
    """
    Create a basic WCS object (with optional shift) and
    append it to the input_dm.meta attribute.

    The WCS object will have a pipeline with the required
    steps to validate against the TweakReg pipeline.

    Parameters
    ----------
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.
    shift_1 : int, optional
        Shift to be applied in the direction of the first axis, by default 0
    shift_2 : int, optional
        Shift to be applied in the direction of the second axis, by default 0
    """

    shape = input_dm.data.shape

    # create necessary transformations
    distortion = Shift(-shift_1) & Shift(-shift_2)
    distortion.bounding_box = ((-0.5, shape[-1] + 0.5), (-0.5, shape[-2] + 0.5))
    tel2sky = pointing.v23tosky(input_dm)

    # create required frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # create pipeline
    pipeline = [
        wcs.Step(detector, distortion),
        wcs.Step(v2v3, tel2sky),
        wcs.Step(world, None),
    ]

    wcs_obj = wcs.WCS(pipeline)

    input_dm.meta["wcs"] = wcs_obj


@pytest.fixture
def base_image():
    """
    Create a base image with a realistic WCSInfo and a WCS.

    Notes
    -----
    The size of the image needs to be relatively large in order for
    the source catalog step to find a reasonable number of sources in the image.

    shift_1 and shift_2 (units in pixel) are used to shift the WCS projection plane.
    """

    def _base_image(shift_1=0, shift_2=0):
        l2 = maker_utils.mk_level2_image(shape=(100, 100))
        l2_im = rdm.ImageModel(l2)
        _create_wcs(l2_im)
        l2_im.meta.wcsinfo.vparity = -1
        return l2_im

    return _base_image
