import numpy as np
import pytest
from astropy.time import Time
from roman_datamodels import stnode
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

    wfi_image = stnode.WfiImage.create_fake_data(shape=shape)
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = "WFI01"
    wfi_image.meta.instrument.optical_element = "F158"
    wfi_image.meta.exposure.type = exptype
    wfi_image.data = np.ones(shape, dtype=np.float32)
    wfi_image.dq = np.zeros(shape, dtype=np.uint32)
    wfi_image.err = np.zeros(shape, dtype=np.float32)
    wfi_image.var_poisson = np.zeros(shape, dtype=np.float32)
    wfi_image.var_rnoise = np.zeros(shape, dtype=np.float32)
    wfi_image.var_flat = np.zeros(shape, dtype=np.float32)

    wfi_image_model = ImageModel(wfi_image)
    wfi_image_model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    wfi_image_model.meta.cal_logs = stnode.CalLogs.create_fake_data()
    flatref = stnode.FlatRef.create_fake_data()
    meta = flatref.meta
    meta["instrument"]["optical_element"] = "F158"
    meta["instrument"]["detector"] = "WFI01"
    flatref["meta"] = meta
    flatref["data"] = np.ones(shape, dtype=np.float32)
    flatref["dq"] = np.zeros(shape, dtype=np.uint32)
    flatref["err"] = (RNG.uniform(size=shape) * 0.05).astype(np.float32)
    flatref_model = FlatRefModel(flatref)

    result = FlatFieldStep.call(wfi_image_model, override_flat=flatref_model)

    assert (result.data == wfi_image.data).all()
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

    shape = (20, 20)

    wfi_image = stnode.WfiImage.create_fake_data(shape=shape)
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = "WFI01"
    wfi_image.meta.instrument.optical_element = "F158"

    wfi_image.meta.exposure.start_time = Time("2020-01-02T11:11:11.110")
    wfi_image.meta.exposure.end_time = Time("2020-01-02T11:33:11.110")

    wfi_image.meta.exposure.type = exptype
    wfi_image_model = ImageModel(wfi_image)
    wfi_image_model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    wfi_image_model.meta.cal_logs = stnode.CalLogs.create_fake_data()

    step = FlatFieldStep()
    ref_file_path = step.get_reference_file(wfi_image_model, "flat")

    wfi_image_model.meta.exposure.start_time = Time("2021-09-01T00:11:11.110")
    wfi_image_model.meta.exposure.end_time = Time("2021-09-01T00:33:11.110")
    ref_file_path_b = step.get_reference_file(wfi_image_model, "flat")
    assert "/".join(ref_file_path.rsplit("/", 1)[1:]) != "/".join(
        ref_file_path_b.rsplit("/", 1)[1:]
    )


@pytest.mark.parametrize("include_var_flat", (True, False))
def test_skip_var_flat(include_var_flat):
    """Test that we don't populate var_flat if requested."""
    wfi_image_model = ImageModel.create_fake_data(shape=(4088, 4088))
    wfi_image_model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    wfi_image_model.meta.cal_logs = stnode.CalLogs.create_fake_data()
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

    wfi_image = stnode.WfiImage.create_fake_data(shape=shape)
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = "WFI01"
    wfi_image.meta.instrument.optical_element = "F158"

    wfi_image.meta.exposure.start_time = Time("2020-02-01T00:00:00.000")
    wfi_image.meta.exposure.end_time = Time("2020-02-01T00:00:05.000")

    wfi_image.meta.exposure.type = exptype
    wfi_image_model = ImageModel(wfi_image)
    wfi_image_model.meta.cal_step = stnode.L2CalStep.create_fake_data()
    wfi_image_model.meta.cal_logs = stnode.CalLogs.create_fake_data()

    result = FlatFieldStep.call(wfi_image_model)
    assert result.meta.cal_step.flat_field == "SKIPPED"
