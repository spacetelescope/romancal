import pytest
from astropy.time import Time
from roman_datamodels.datamodels import FlatRefModel, WfiImageModel
from stpipe import crds_client

from romancal.stpipe import RomanPipeline, RomanStep


@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_open_model(step_class, tmp_path):
    """
    Test that the class is properly hooked up to datamodels.open.
    More comprehensive tests can be found in romancal.datamodels.tests,
    this is just a smoke test of the integration.
    """
    file_path = tmp_path / "test.asdf"

    WfiImageModel.make_default(shape=(8, 20, 20), filepath=file_path)

    step = step_class()
    with step.open_model(file_path) as model:
        assert model.meta.telescope == "ROMAN"


@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_get_reference_file(step_class):
    """
    Test that CRDS is properly integrated.
    """
    model = WfiImageModel.make_default(shape=(8, 20, 20))
    # This will be brittle while we're using the dev server.
    # If this test starts failing mysteriously, check the
    # metadata values against the flat rmap.
    model.meta.instrument.optical_element = "F158"
    model.meta.exposure.start_time = Time("2021-01-01T12:00:00")

    step = step_class()
    reference_path = step.get_reference_file(model, "flat")

    with step.open_model(reference_path) as reference_model:
        assert isinstance(reference_model, FlatRefModel)


@pytest.mark.skip(reason="There are no grism flats.")
@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_get_reference_file_spectral(step_class):
    """
    Test that CRDS is properly integrated.
    """
    model = WfiImageModel.make_default(shape=(8, 20, 20))
    # This will be brittle while we're using the dev server.
    # If this test starts failing mysteriously, check the
    # metadata values against the flat rmap.
    model.meta.instrument.optical_element = "GRISM"
    model.meta.exposure.start_time = Time("2021-01-01T12:00:00")

    step = step_class()
    reference_path = step.get_reference_file(model, "flat")

    with step.open_model(reference_path) as reference_model:
        assert isinstance(reference_model, FlatRefModel)
        assert reference_model.meta.instrument.optical_element == "GRISM"


def test_log_messages(tmp_path):
    class LoggingStep(RomanStep):
        def process(self):
            self.log.warning("Splines failed to reticulate")
            return WfiImageModel.make_default(shape=(8, 20, 20))

    result = LoggingStep().run()
    assert any("Splines failed to reticulate" in l for l in result.cal_logs)


def test_crds_meta():
    """Test that context and software versions are set"""

    class ReflectStep(RomanStep):
        def process(self, input):
            return input

    im = WfiImageModel.make_default(shape=(8, 20, 20))
    im.meta.ref_file.crds.sw_version = "junkversion"
    im.meta.ref_file.crds.context_used = "junkcontext"

    result = ReflectStep.call(im)

    assert result.meta.ref_file.crds.sw_version == crds_client.get_svn_version()
    assert result.meta.ref_file.crds.context_used == crds_client.get_context_used(
        result.crds_observatory
    )
