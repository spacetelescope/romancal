"""Fixtures for linearity tests."""

import numpy as np
import pytest
from astropy.time import Time
from roman_datamodels.datamodels import ScienceRawModel

from romancal.dq_init import DQInitStep


@pytest.fixture(scope="function")
def setup_ramp_for_linearity():
    """Set up ramp model ready for linearity correction.

    Returns a factory function that creates a ramp model with metadata
    populated, DQ initialized, and read_pattern set.
    """

    def _model(shape):
        model = ScienceRawModel.create_fake_data(shape=shape)
        model.meta.exposure.start_time = Time(
            "2024-01-03T00:00:00.0", format="isot", scale="utc"
        )
        model.meta.instrument.name = "WFI"
        model.meta.instrument.detector = "WFI01"
        model.meta.instrument.optical_element = "F158"
        model.meta["guide_star"]["window_xstart"] = 1012
        model.meta["guide_star"]["window_xsize"] = 16
        model.meta.exposure.type = "WFI_IMAGE"
        model.data = np.zeros(shape, dtype=np.uint16)
        model.meta.exposure.read_pattern = [[i + 1] for i in range(shape[0])]

        result = DQInitStep.call(model)
        return result

    return _model
