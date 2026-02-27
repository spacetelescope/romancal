import json
import logging

import pytest
from astropy.time import Time
from roman_datamodels.datamodels import FlatRefModel, ImageModel
from stpipe import crds_client

import romancal
from romancal.datamodels import ModelLibrary
from romancal.flatfield import FlatFieldStep
from romancal.stpipe import RomanPipeline, RomanStep


@pytest.mark.parametrize("is_container", [True, False])
@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_open_model(step_class, tmp_path, is_container, base_image):
    """
    Test that the class is properly hooked up to datamodels.open.
    More comprehensive tests can be found in romancal.datamodels.tests,
    this is just a smoke test of the integration.
    """
    file_path = tmp_path / "test.asdf"

    im = base_image()
    im.save(file_path)

    if is_container:
        asn = {
            "asn_pool": "none",
            "products": [
                {
                    "members": [
                        {
                            "exptype": "science",
                            "expname": "test.asdf",
                        }
                    ],
                },
            ],
        }
        asn_path = tmp_path / "test.json"
        with open(asn_path, "w") as f:
            json.dump(asn, f)
        test_file_path = asn_path
    else:
        test_file_path = file_path

    step = step_class()
    with step.open_model(test_file_path) as model:
        if is_container:
            assert isinstance(model, ModelLibrary)
            assert model.crds_observatory == "roman"
            assert model.get_crds_parameters() is not None
        else:
            assert isinstance(model, ImageModel)
            assert model.meta.telescope == "ROMAN"


@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_get_reference_file(step_class, base_image):
    """
    Test that CRDS is properly integrated.
    """
    im = base_image()
    # This will be brittle while we're using the dev server.
    # If this test starts failing mysteriously, check the
    # metadata values against the flat rmap.
    im.meta.instrument.optical_element = "F158"
    im.meta.exposure.start_time = Time("2024-01-01T12:00:00")

    step = step_class()
    reference_path = step.get_reference_file(im, "flat")

    with step.open_model(reference_path) as reference_model:
        assert isinstance(reference_model, FlatRefModel)


@pytest.mark.skip(reason="There are no grism flats.")
@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_get_reference_file_spectral(step_class, base_image):
    """
    Test that CRDS is properly integrated.
    """
    im = base_image()
    # This will be brittle while we're using the dev server.
    # If this test starts failing mysteriously, check the
    # metadata values against the flat rmap.
    im.meta.instrument.optical_element = "GRISM"
    im.meta.exposure.start_time = Time("2024-01-01T12:00:00")

    step = step_class()
    reference_path = step.get_reference_file(im, "flat")

    with step.open_model(reference_path) as reference_model:
        assert isinstance(reference_model, FlatRefModel)
        assert reference_model.meta.instrument.optical_element == "GRISM"


def test_log_messages(tmp_path, base_image):
    LOGGER_NAME = "test_log_messages"
    logger = logging.getLogger(LOGGER_NAME)

    class LoggingStep(RomanStep):
        def process(self):
            logger.warning("Splines failed to reticulate")
            return base_image()

        @staticmethod
        def get_stpipe_loggers():
            return (*RomanStep.get_stpipe_loggers(), LOGGER_NAME)

    result = LoggingStep().run()
    assert any("Splines failed to reticulate" in l for l in result.meta.cal_logs)


def test_crds_meta(base_image):
    """Test that context and software versions are set"""

    im = base_image()
    result = FlatFieldStep.call(im)

    assert result.meta.ref_file.crds.version == crds_client.get_svn_version()
    assert result.meta.ref_file.crds.context == crds_client.get_context_used(
        result.crds_observatory
    )


def test_calibration_software_version(base_image):
    """Test that calibration_software_version is updated when a step is run"""

    class NullStep(RomanStep):
        def process(self, input):
            return input

    im = base_image()
    im.meta.calibration_software_version = "junkversion"

    result = NullStep.call(im)

    assert result.meta.calibration_software_version == romancal.__version__
