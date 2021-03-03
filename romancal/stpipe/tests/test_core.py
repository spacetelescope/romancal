import os

import asdf
import pytest

from romancal.datamodels import ImageModel, FlatModel
from romancal.stpipe import RomanPipeline, RomanStep


@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_open_model(step_class, tmp_path):
    """
    Test that the class is properly hooked up to datamodels.open.
    More comprehensive tests can be found in romancal.datamodels.tests,
    this is just a smoke test of the integration.
    """
    file_path = tmp_path / "test.asdf"

    with asdf.AsdfFile() as af:
        af["meta"] = {"telescope": "a dashing monocle"}
        af.write_to(file_path)

    step = step_class()
    with step.open_model(file_path) as model:
        assert model.meta.telescope == "a dashing monocle"


@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently available outside the internal network"
)
@pytest.mark.parametrize("step_class", [RomanPipeline, RomanStep])
def test_get_reference_file(step_class):
    """
    Test that CRDS is properly integrated.
    """
    model = ImageModel()
    # This will be brittle while we're using the dev server.
    # If this test starts failing mysteriously, check the
    # metadata values against the flat rmap.
    model.meta.instrument.name = "WFI"
    model.meta.instrument.detector = "WFI01"
    model.meta.instrument.optical_element = "F158"
    model.meta.observation.date = "2020-01-01"
    model.meta.observation.time = "00:00:00"

    step = step_class()
    reference_path = step.get_reference_file(model, "flat")

    with step.open_model(reference_path) as reference_model:
        assert isinstance(reference_model, FlatModel)
        assert reference_model.meta.instrument.name == 'WFI'
        assert reference_model.meta.instrument.detector == 'WFI01'
        assert reference_model.meta.instrument.optical_element == 'F158'
        assert reference_model.data.shape == (4096, 4096)
