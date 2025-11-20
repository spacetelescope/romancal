import numpy as np
import pytest
from astropy.time import Time
from roman_datamodels.datamodels import FlatRefModel, ImageModel

from romancal.flatfield import FlatFieldStep

RNG = np.random.default_rng(172)


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_flatfield_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a FLAT reffile"""

    shape = (20, 20)

    wfi_image_model = ImageModel.create_fake_data(shape=shape)
    wfi_image_model.meta.exposure.start_time = Time(
        "2024-01-03T00:00:00.0", format="isot", scale="utc"
    )
    wfi_image_model.meta.instrument.name = instrument
    wfi_image_model.meta.instrument.detector = "WFI01"
    wfi_image_model.meta.instrument.optical_element = "F158"
    wfi_image_model.meta.exposure.type = exptype
    wfi_image_model.data = np.ones(shape, dtype=np.float32)
    wfi_image_model.dq = np.zeros(shape, dtype=np.uint32)
    wfi_image_model.err = np.zeros(shape, dtype=np.float32)
    wfi_image_model.var_poisson = np.zeros(shape, dtype=np.float32)
    wfi_image_model.var_rnoise = np.zeros(shape, dtype=np.float32)
    wfi_image_model.var_flat = np.zeros(shape, dtype=np.float32)
    wfi_image_model.meta.cal_step = {}
    for step_name in wfi_image_model.schema_info("required")["roman"]["meta"][
        "cal_step"
    ]["required"].info:
        wfi_image_model.meta.cal_step[step_name] = "INCOMPLETE"
    wfi_image_model.meta.cal_logs = []

    flatref_model = FlatRefModel.create_fake_data()
    flatref_model["meta"]["instrument"]["optical_element"] = "F158"
    flatref_model["meta"]["instrument"]["detector"] = "WFI01"
    flatref_model["data"] = np.ones(shape, dtype=np.float32)
    flatref_model["dq"] = np.zeros(shape, dtype=np.uint32)
    flatref_model["err"] = (RNG.uniform(size=shape) * 0.05).astype(np.float32)

    result = FlatFieldStep.call(wfi_image_model, override_flat=flatref_model)

    assert (result.data == wfi_image_model.data).all()
    assert result.var_flat.shape == shape
    assert result.meta.cal_step.flat_field == "COMPLETE"

    # test that the step is skipped if the reference file is N/A
    result = FlatFieldStep.call(wfi_image_model, override_flat="N/A")

    assert result.meta.cal_step.flat_field == "SKIPPED"


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_crds_temporal_match(instrument, exptype):
    """Test that the basic inferface works for data requiring a FLAT reffile"""

    from crds.core.exceptions import CrdsLookupError

    shape = (20, 20)

    wfi_image_model = ImageModel.create_fake_data(shape=shape)
    wfi_image_model.meta.instrument.name = instrument
    wfi_image_model.meta.instrument.detector = "WFI01"
    wfi_image_model.meta.instrument.optical_element = "F158"

    wfi_image_model.meta.exposure.start_time = Time("2025-01-02T11:11:11.110")
    wfi_image_model.meta.exposure.end_time = Time("2025-01-02T11:33:11.110")

    wfi_image_model.meta.exposure.type = exptype
    wfi_image_model.meta.cal_step = {}
    for step_name in wfi_image_model.schema_info("required")["roman"]["meta"][
        "cal_step"
    ]["required"].info:
        wfi_image_model.meta.cal_step[step_name] = "INCOMPLETE"

    wfi_image_model.meta.cal_logs = []

    step = FlatFieldStep()
    ref_file_path = step.get_reference_file(wfi_image_model, "flat")
    assert 'roman_wfi_flat' in ref_file_path

    # As of 20251120, there are no references with multiple USEAFTER.
    # However, picking a date earlier than the first USEAFTER
    # and expect that the retrieval will fail.
    wfi_image_model.meta.exposure.start_time = Time("2021-09-01T00:11:11.110")
    wfi_image_model.meta.exposure.end_time = Time("2021-09-01T00:33:11.110")
    with pytest.raises(CrdsLookupError):
        _ = step.get_reference_file(wfi_image_model, "flat")



@pytest.mark.parametrize("include_var_flat", (True, False))
def test_skip_var_flat(include_var_flat):
    """Test that we don't populate var_flat if requested."""
    wfi_image_model = ImageModel.create_fake_data(shape=(4088, 4088))
    wfi_image_model.meta.exposure.start_time = Time(
        "2024-01-03T00:00:00.0", format="isot", scale="utc"
    )
    wfi_image_model.meta.cal_step = {}
    for step_name in wfi_image_model.schema_info("required")["roman"]["meta"][
        "cal_step"
    ]["required"].info:
        wfi_image_model.meta.cal_step[step_name] = "INCOMPLETE"
    wfi_image_model.meta.cal_logs = []
    result = FlatFieldStep.call(wfi_image_model, include_var_flat=include_var_flat)
    assert hasattr(result, "var_flat") == include_var_flat


@pytest.mark.parametrize(
    "instrument",
    [
        "WFI",
    ],
)
@pytest.mark.parametrize(
    "exptype",
    [
        "WFI_GRISM",
        "WFI_PRISM",
    ],
)
def test_spectroscopic_skip(instrument, exptype):
    shape = (20, 20)

    wfi_image_model = ImageModel.create_fake_data(shape=shape)
    wfi_image_model.meta.instrument.name = instrument
    wfi_image_model.meta.instrument.detector = "WFI01"
    wfi_image_model.meta.instrument.optical_element = "F158"

    wfi_image_model.meta.exposure.start_time = Time("2024-02-01T00:00:00.000")
    wfi_image_model.meta.exposure.end_time = Time("2024-02-01T00:00:05.000")

    wfi_image_model.meta.exposure.type = exptype
    wfi_image_model.meta.cal_step = {}
    for step_name in wfi_image_model.schema_info("required")["roman"]["meta"][
        "cal_step"
    ]["required"].info:
        wfi_image_model.meta.cal_step[step_name] = "INCOMPLETE"
    wfi_image_model.meta.cal_logs = []

    result = FlatFieldStep.call(wfi_image_model)
    assert result.meta.cal_step.flat_field == "SKIPPED"
