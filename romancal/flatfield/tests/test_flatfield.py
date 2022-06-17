import os

import pytest
import numpy as np

from roman_datamodels import stnode
from roman_datamodels.datamodels import ImageModel, FlatRefModel
from roman_datamodels.testing import utils as testutil
from romancal.flatfield import FlatFieldStep
from astropy.time import Time


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ]
)
@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently \
    available outside the internal network"
)
def test_flatfield_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a FLAT reffile"""

    shape = (20, 20)

    wfi_image = testutil.mk_level2_image(shape=shape)
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = 'WFI01'
    wfi_image.meta.instrument.optical_element = 'F158'
    wfi_image.meta.exposure.type = exptype
    wfi_image.data = np.ones(shape, dtype=np.float32)
    wfi_image.dq = np.zeros(shape, dtype=np.uint32)
    wfi_image.err = np.zeros(shape, dtype=np.float32)
    wfi_image.var_poisson = np.zeros(shape, dtype=np.float32)
    wfi_image.var_rnoise = np.zeros(shape, dtype=np.float32)
    wfi_image.var_flat = np.zeros(shape, dtype=np.float32)
    wfi_image_model = ImageModel(wfi_image)
    flatref = stnode.FlatRef()
    meta = {}
    testutil.add_ref_common(meta)
    meta['instrument']['optical_element'] = 'F158'
    meta['instrument']['detector'] = 'WFI01'
    meta['reftype'] = 'FLAT'
    flatref['meta'] = meta
    flatref['data'] = np.ones(shape, dtype=np.float32)
    flatref['dq'] = np.zeros(shape, dtype=np.uint16)
    flatref['err'] = (np.random.random(shape) * 0.05).astype(np.float32)
    flatref_model = FlatRefModel(flatref)

    result = FlatFieldStep.call(wfi_image_model, override_flat=flatref_model)

    assert (result.data == wfi_image.data).all()
    assert result.var_flat.shape == shape
    assert result.meta.cal_step.flat_field == 'COMPLETE'


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ]
)
@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently \
    available outside the internal network"
)
def test_crds_temporal_match(instrument, exptype):
    """Test that the basic inferface works for data requiring a FLAT reffile"""

    wfi_image = testutil.mk_level2_image()
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = 'WFI01'
    wfi_image.meta.instrument.optical_element = 'F158'

    wfi_image.meta.exposure.start_time = Time('2020-01-02T11:11:11.110')
    wfi_image.meta.exposure.end_time = Time('2020-01-02T11:33:11.110')

    wfi_image.meta.exposure.type = exptype
    wfi_image_model = ImageModel(wfi_image)

    step = FlatFieldStep()
    ref_file_path = step.get_reference_file(wfi_image_model, "flat")

    wfi_image_model.meta.exposure.start_time = Time('2021-09-01T00:11:11.110')
    wfi_image_model.meta.exposure.end_time = Time('2021-09-01T00:33:11.110')
    ref_file_path_b = step.get_reference_file(wfi_image_model, "flat")
    assert ("/".join(ref_file_path.rsplit("/", 1)[1:])) != \
           ("/".join(ref_file_path_b.rsplit("/", 1)[1:]))


@pytest.mark.parametrize(
    "instrument", ["WFI", ]
)
@pytest.mark.parametrize(
    "exptype", ["WFI_GRISM", "WFI_PRISM", ]
)
@pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Roman CRDS servers are not currently \
    available outside the internal network"
)
# Test that spectroscopic exposure types will skip flat field step
def test_spectroscopic_skip(instrument, exptype):
    wfi_image = testutil.mk_level2_image()
    wfi_image.meta.instrument.name = instrument
    wfi_image.meta.instrument.detector = 'WFI01'
    wfi_image.meta.instrument.optical_element = 'F158'

    wfi_image.meta.exposure.start_time = Time('2020-02-01T00:00:00.000')
    wfi_image.meta.exposure.end_time = Time('2020-02-01T00:00:05.000')

    wfi_image.meta.exposure.type = exptype
    wfi_image_model = ImageModel(wfi_image)

    result = FlatFieldStep.call(wfi_image_model)
    assert result.meta.cal_step.flat_field == 'SKIPPED'
