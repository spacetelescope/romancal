import json
import logging

import numpy as np
import pytest
from astropy.time import Time
from roman_datamodels.datamodels import FlatRefModel, ImageModel, MosaicModel
from roman_datamodels.dqflags import pixel
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


class MockStepClass(RomanStep):
    """Minimal subclass to test RomanStep methods."""

    def process(self, input):
        return input


@pytest.mark.parametrize(
    "model_type, model_class",
    [("imagemodel", ImageModel), ("mosaicmodel", MosaicModel)],
)
@pytest.mark.parametrize(
    "data_val, dq_val, expected_frac",
    [
        # all pixels good (dq=0)
        (10.0, 0, 1.0),
        # all pixels flagged do_not_use (fraction should be 0.0)
        (10.0, pixel.DO_NOT_USE, 0.0),
        # half pixels flagged do_not_use (0.5)
        (5.0, "half_bad", 0.5),
        # other flags set (saturated), but not do_not_use (fraction should be 1.0)
        (10.0, pixel.SATURATED | pixel.HOT, 1.0),
    ],
)
def test_populate_statistics(model_type, model_class, data_val, dq_val, expected_frac):
    """Test statistics population."""
    step = MockStepClass()
    shape = (10, 10)

    model = model_class.create_minimal()

    model.data = np.full(shape, data_val, dtype=np.float32)

    if dq_val == "half_bad":
        dq = np.zeros(shape, dtype=np.uint32)
        dq.view().flat[0:50] = pixel.DO_NOT_USE
        model.dq = dq
    else:
        model.dq = np.full(shape, dq_val, dtype=np.uint32)

    # inject a nan to verify nan-resistance (nanmedian and mad_std)
    model.data[0, 0] = np.nan

    step.populate_statistics(model)

    assert model.meta.statistics.image_median == data_val

    assert model.meta.statistics.image_rms == 0.0

    assert model.meta.statistics.good_pixel_fraction == pytest.approx(
        expected_frac, abs=1e-6
    )

    # check the placeholder logic (-1.0)
    assert model.meta.statistics.zodiacal_light == -1.0


def test_statistics_handled_in_finalize():
    """Verify finalize_result triggers the statistics population."""
    step = MockStepClass()
    shape = (10, 10)

    model = ImageModel.create_minimal()

    model.data = np.ones(shape, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)

    # finalize_result should trigger the populate_statistics call
    step.finalize_result(model, [])

    assert hasattr(model.meta, "statistics")
    assert model.meta.statistics.image_median == 1.0
    assert model.meta.statistics.good_pixel_fraction == 1.0


def test_statistics_graceful_exit_no_data():
    """Ensure we don't crash if data is None."""
    step = MockStepClass()

    model = ImageModel.create_minimal()

    model.data = None

    step.populate_statistics(model)

    # check that meta.statistics wasn't created
    assert not hasattr(model.meta, "statistics")
