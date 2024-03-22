import numpy as np
import pytest
from astropy import units as u
from roman_datamodels import maker_utils, stnode
from roman_datamodels.datamodels import MaskRefModel, ScienceRawModel
from roman_datamodels.dqflags import pixel

from romancal.dq_init import DQInitStep
from romancal.dq_init.dq_initialization import do_dqinit

RNG = np.random.default_rng(83)

# Set parameters for multiple runs of data
args = "xstart, ystart, xsize, ysize, ngroups, instrument, exp_type"
test_data = [(1, 1, 2048, 2048, 2, "WFI", "WFI_IMAGE")]
ids = ["RampModel"]


@pytest.mark.parametrize(args, test_data, ids=ids)
def test_dq_im(xstart, ystart, xsize, ysize, ngroups, instrument, exp_type):
    """Check that PIXELDQ is initialized with the information from the reference file.
    test that a flagged value in the reference file flags the PIXELDQ array"""

    csize = (ngroups, ysize, xsize)

    # create raw input data for step
    dm_ramp = maker_utils.mk_ramp(shape=csize)
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    dq = np.zeros(csize[1:], dtype=np.uint32)

    # edit reference file with known bad pixel values
    dq[100, 100] = 2  # Saturated pixel
    dq[200, 100] = 4  # Jump detected pixel
    dq[300, 100] = 8  # Dropout
    dq[400, 100] = 32  # Persistence
    dq[500, 100] = 1  # Do_not_use
    dq[600, 100] = 16  # guide window affected data
    dq[100, 200] = 3  # Saturated pixel + do not use
    dq[200, 200] = 5  # Jump detected pixel + do not use
    dq[300, 200] = 9  # Dropout + do not use
    dq[400, 200] = 33  # Persistence + do not use

    # write mask model
    ref_data = maker_utils.mk_mask(shape=csize[1:])

    # Copy in maskmodel elemnts
    ref_data["dq"] = dq

    ref_data["meta"]["instrument"]["name"] = instrument

    # run do_dqinit
    outfile = do_dqinit(dm_ramp, ref_data)
    dqdata = outfile["pixeldq"]

    # assert that the pixels read back in match the mapping from ref data to
    # science data
    assert dqdata[100, 100] == pixel.SATURATED
    assert dqdata[200, 100] == pixel.JUMP_DET
    assert dqdata[300, 100] == pixel.DROPOUT
    assert dqdata[400, 100] == pixel.PERSISTENCE
    assert dqdata[500, 100] == pixel.DO_NOT_USE
    assert dqdata[600, 100] == pixel.GW_AFFECTED_DATA
    assert dqdata[100, 200] == pixel.SATURATED + pixel.DO_NOT_USE
    assert dqdata[200, 200] == pixel.JUMP_DET + pixel.DO_NOT_USE
    assert dqdata[300, 200] == pixel.DROPOUT + pixel.DO_NOT_USE
    assert dqdata[400, 200] == pixel.PERSISTENCE + pixel.DO_NOT_USE


def test_groupdq():
    """
    Check that GROUPDQ extension is added to the data and all values are
    initialized to zero.
    """

    # size of integration
    instrument = "WFI"
    ngroups = 5
    xsize = 1032
    ysize = 1024
    csize = (ngroups, ysize, xsize)

    # create raw input data for step
    dm_ramp = maker_utils.mk_ramp(shape=csize)
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    ref_data = maker_utils.mk_mask(shape=csize[1:])
    ref_data["meta"]["instrument"]["name"] = instrument

    # run the correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that GROUPDQ was created and initialized to zero
    groupdq = outfile.groupdq

    np.testing.assert_array_equal(
        np.full((ngroups, ysize, xsize), 0, dtype=int),
        groupdq,
        err_msg="groupdq not initialized to zero",
    )


def test_err():
    """Check that a 3-D ERR array is initialized and all values are zero."""

    # size of integration
    instrument = "WFI"
    ngroups = 5
    xsize = 1032
    ysize = 1024
    csize = (ngroups, ysize, xsize)

    # create raw input data for step
    dm_ramp = maker_utils.mk_ramp(shape=(ngroups, ysize, xsize))
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    ref_data = maker_utils.mk_mask(shape=csize[1:])
    ref_data["meta"]["instrument"]["name"] = instrument

    # run correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that ERR array was created and initialized to zero
    errarr = outfile.err

    assert errarr.ndim == 3  # check that output err array is 3-D
    assert np.all(errarr == 0)  # check that values are 0


def test_dq_add1_groupdq():
    """
    Test if the dq_init code set the groupdq flag on the first
    group to 'do_not_use' by adding 1 to the flag, not overwriting to 1
    Also test whether two flags on the same pixel are added together.
    """

    # size of integration
    instrument = "WFI"
    ngroups = 5
    xsize = 1032
    ysize = 1024
    csize = (ngroups, ysize, xsize)

    # create raw input data for step
    dm_ramp = maker_utils.mk_ramp(shape=(ngroups, ysize, xsize))
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    dq = np.zeros(csize[1:], dtype=np.uint32)

    # write reference file with known bad pixel values
    dq[505, 505] = 1  # Do_not_use
    dq[400, 500] = 3  # do_not_use and saturated pixel

    # write mask model
    ref_data = maker_utils.mk_mask(shape=csize[1:])
    ref_data["meta"]["instrument"]["name"] = instrument

    # Copy in maskmodel elemnts
    ref_data["dq"] = dq

    # set a flag in the pixel dq
    dm_ramp.pixeldq[505, 505] = 4  # Jump detected pixel

    # run correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # test if pixels in pixeldq were incremented in value by 1
    # check that previous dq flag is added to mask value
    assert outfile.pixeldq[505, 505] == pixel.JUMP_DET + pixel.DO_NOT_USE
    # check two flags propagate correctly
    assert outfile.pixeldq[400, 500] == pixel.SATURATED + pixel.DO_NOT_USE


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dqinit_step_interface(instrument, exptype):
    """Test that the basic inferface works for data requiring a DQ reffile"""

    # Set test size
    shape = (2, 20, 20)

    # include an extra item not defined in the schema
    extra_key = "foo_extra"
    extra_value = [1, 2, 3]

    # Create test science raw model
    wfi_sci_raw = maker_utils.mk_level1_science_raw(shape=shape)
    wfi_sci_raw.meta.instrument.name = instrument
    wfi_sci_raw.meta.instrument.detector = "WFI01"
    wfi_sci_raw.meta.instrument.optical_element = "F158"
    wfi_sci_raw.meta["guidestar"]["gw_window_xstart"] = 1012
    wfi_sci_raw.meta["guidestar"]["gw_window_xsize"] = 16
    wfi_sci_raw.meta.exposure.type = exptype
    wfi_sci_raw.data = u.Quantity(
        np.ones(shape, dtype=np.uint16), u.DN, dtype=np.uint16
    )
    wfi_sci_raw[extra_key] = extra_value
    wfi_sci_raw_model = ScienceRawModel(wfi_sci_raw)

    # Create mask model
    maskref = stnode.MaskRef()
    meta = maker_utils.mk_ref_common("MASK")
    meta["instrument"]["optical_element"] = "F158"
    meta["instrument"]["detector"] = "WFI01"
    maskref["meta"] = meta
    maskref["data"] = np.ones(shape[1:], dtype=np.float32)
    maskref["dq"] = np.zeros(shape[1:], dtype=np.uint32)
    maskref["err"] = (RNG.uniform(size=shape[1:]) * 0.05).astype(np.float32)
    maskref_model = MaskRefModel(maskref)

    # Perform Data Quality application step
    result = DQInitStep.call(wfi_sci_raw_model, override_mask=maskref_model)

    # Test dq_init results
    assert (result.data == wfi_sci_raw.data).all()
    assert result.pixeldq.shape == shape[1:]
    assert result.meta.cal_step.dq_init == "COMPLETE"
    assert result.data.dtype == np.float32
    assert result.err.dtype == np.float32
    assert result.pixeldq.dtype == np.uint32
    assert result.groupdq.dtype == np.uint8
    # check that extra value came through
    assert result[extra_key] == extra_value


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dqinit_refpix(instrument, exptype):
    """Test that the basic inferface works for data requiring a DQ reffile"""

    # Set test size
    shape = (2, 20, 20)

    # Create test science raw model
    wfi_sci_raw = maker_utils.mk_level1_science_raw(shape=shape)
    wfi_sci_raw.meta.instrument.name = instrument
    wfi_sci_raw.meta.instrument.detector = "WFI01"
    wfi_sci_raw.meta.instrument.optical_element = "F158"
    wfi_sci_raw.meta["guidestar"]["gw_window_xstart"] = 1012
    wfi_sci_raw.meta["guidestar"]["gw_window_xsize"] = 16
    wfi_sci_raw.meta.exposure.type = exptype
    wfi_sci_raw.data = u.Quantity(
        np.ones(shape, dtype=np.uint16), u.DN, dtype=np.uint16
    )
    wfi_sci_raw_model = ScienceRawModel(wfi_sci_raw)

    # Create mask model
    maskref = stnode.MaskRef()
    meta = maker_utils.mk_ref_common("MASK")
    meta["instrument"]["optical_element"] = "F158"
    meta["instrument"]["detector"] = "WFI01"
    maskref["meta"] = meta
    maskref["data"] = np.ones(shape[1:], dtype=np.float32)
    maskref["dq"] = np.zeros(shape[1:], dtype=np.uint32)
    maskref["err"] = (RNG.uniform(size=shape[1:]) * 0.05).astype(np.float32)
    maskref_model = MaskRefModel(maskref)

    # Perform Data Quality application step
    result = DQInitStep.call(wfi_sci_raw_model, override_mask=maskref_model)

    # check if reference pixels are correct
    assert result.data.shape == (2, 20, 20)  # no pixels should be trimmed
    assert result.amp33.value.shape == (2, 4096, 128)
    assert result.border_ref_pix_right.shape == (2, 20, 4)
    assert result.border_ref_pix_left.shape == (2, 20, 4)
    assert result.border_ref_pix_top.shape == (2, 4, 20)
    assert result.border_ref_pix_bottom.shape == (2, 4, 20)
    assert result.dq_border_ref_pix_right.shape == (20, 4)
    assert result.dq_border_ref_pix_left.shape == (20, 4)
    assert result.dq_border_ref_pix_top.shape == (4, 20)
    assert result.dq_border_ref_pix_bottom.shape == (4, 20)


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dqinit_resultantdq(instrument, exptype):
    """Test that the basic inferface works for data requiring a DQ reffile"""

    # Set test size
    shape = (2, 20, 20)

    # Create test science raw model
    wfi_sci_raw = maker_utils.mk_level1_science_raw(shape=shape, dq=True)
    wfi_sci_raw.meta.instrument.name = instrument
    wfi_sci_raw.meta.instrument.detector = "WFI01"
    wfi_sci_raw.meta.instrument.optical_element = "F158"
    wfi_sci_raw.meta["guidestar"]["gw_window_xstart"] = 1012
    wfi_sci_raw.meta["guidestar"]["gw_window_xsize"] = 16
    wfi_sci_raw.meta.exposure.type = exptype
    wfi_sci_raw.resultantdq[1, 12, 12] = pixel["DROPOUT"]
    wfi_sci_raw.data = u.Quantity(
        np.ones(shape, dtype=np.uint16), u.DN, dtype=np.uint16
    )
    wfi_sci_raw_model = ScienceRawModel(wfi_sci_raw)

    # Create mask model
    maskref = stnode.MaskRef()
    meta = maker_utils.mk_ref_common("MASK")
    meta["instrument"]["optical_element"] = "F158"
    meta["instrument"]["detector"] = "WFI01"
    maskref["meta"] = meta
    maskref["data"] = np.ones(shape[1:], dtype=np.float32)
    maskref["dq"] = np.zeros(shape[1:], dtype=np.uint32)
    maskref["err"] = (RNG.uniform(size=shape[1:]) * 0.05).astype(np.float32)
    maskref_model = MaskRefModel(maskref)

    # Perform Data Quality application step
    result = DQInitStep.call(wfi_sci_raw_model, override_mask=maskref_model)

    # check that the resultantdq is present in the raw model
    assert hasattr(wfi_sci_raw_model, "resultantdq")
    # check that the resultantdq is not copied to the result
    assert not hasattr(result, "resultantdq")
    # check to see the resultantdq is the correct shape
    assert wfi_sci_raw_model.resultantdq.shape == shape
    # check to see the resultantdq & groupdq have the correct value
    assert wfi_sci_raw_model.resultantdq[1, 12, 12] == pixel["DROPOUT"]
    assert result.groupdq[1, 12, 12] == pixel["DROPOUT"]


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_dqinit_getbestref(instrument, exptype):
    """Test that the step is skipped if the CRDS reffile is N/A"""

    # Set test size
    shape = (2, 20, 20)

    # Create test science raw model
    wfi_sci_raw = maker_utils.mk_level1_science_raw(shape=shape)
    wfi_sci_raw.meta.instrument.name = instrument
    wfi_sci_raw.meta.instrument.detector = "WFI01"
    wfi_sci_raw.meta.instrument.optical_element = "F158"
    wfi_sci_raw.meta["guidestar"]["gw_window_xstart"] = 1012
    wfi_sci_raw.meta["guidestar"]["gw_window_xsize"] = 16
    wfi_sci_raw.meta.exposure.type = exptype
    wfi_sci_raw.data = u.Quantity(
        np.ones(shape, dtype=np.uint16), u.DN, dtype=np.uint16
    )
    wfi_sci_raw_model = ScienceRawModel(wfi_sci_raw)

    # Perform Data Quality application step
    result = DQInitStep.call(wfi_sci_raw_model, override_mask="N/A")
    assert result.meta.cal_step.dq_init == "SKIPPED"
