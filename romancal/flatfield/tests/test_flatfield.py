import pytest
import numpy as np

from roman_datamodels import util, stnode
from roman_datamodels.datamodels import ImageModel, FlatRefModel
from roman_datamodels import test_utils as testutil
from romancal.flatfield import FlatFieldStep


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ]
)
@pytest.mark.skip(reason="modifying reference_file_types caused other tests to fail")
def test_flatfield_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a FLAT reffile"""

    shape = (20, 20)

    wfi_image = testutil.mk_level2_image(no_arrays=True)
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
    wfi_image.area = np.ones(shape, dtype=np.float32)
    wfi_image_model = ImageModel(wfi_image)
    #wfi_image_model.data[0, 0] = np.nan
    # data = datamodels.ImageModel(shape)
    # data.meta.instrument.name = instrument
    # data.meta.exposure.type = exptype
    # data.meta.subarray.xstart = 1
    # data.meta.subarray.ystart = 1
    # data.meta.subarray.xsize = shape[1]
    # data.meta.subarray.ysize = shape[0]

    flatref = stnode.FlatRef()
    meta = {}
    testutil.add_ref_common(meta)
    meta['instrument']['optical_element'] = 'F158'
    meta['instrument']['dectector'] = 'WFI01'
    meta['reftype'] = 'FLAT'
    flatref['meta'] = meta
    flatref['data'] = np.ones(shape, dtype=np.float32)
    flatref['dq'] = np.zeros(shape, dtype=np.uint16)
    #flatref.data[0, 0] = np.nan
    flatref['err'] = (np.random.random(shape) * 0.05).astype(np.float32)
    flatref_model = FlatRefModel(flatref)
    # flat = datamodels.FlatModel(shape)
    # flat.meta.instrument.name = instrument
    # flat.meta.subarray.xstart = 1
    # flat.meta.subarray.ystart = 1
    # flat.meta.subarray.xsize = shape[1]
    # flat.meta.subarray.ysize = shape[0]
    # flat.data += 1
    # flat.data[0,0] = np.nan
    # flat.err = np.random.random(shape) * 0.05

    # override class attribute so only the `flat` type needs to be overriden
    # in the step call.
    # FlatFieldStep.reference_file_types = ["flat"]
    result = FlatFieldStep.call(wfi_image_model, override_flat=flatref_model)

    assert (result.data == wfi_image.data).all()
    assert result.var_flat.shape == shape
    assert result.meta.cal_step.flat_field == 'COMPLETE'
