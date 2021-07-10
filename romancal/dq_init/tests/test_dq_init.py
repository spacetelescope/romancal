from stdatamodels.validate import ValidationWarning
import numpy as np
import pytest
import warnings

from romancal.dq_init import DQInitStep
from romancal.dq_init.dq_initialization import do_dqinit
from romancal.stpipe import RomanStep

from roman_datamodels import stnode, table_definitions, dqflags
from roman_datamodels.testing.factories import _random_dq_def
from roman_datamodels.datamodels import MaskRefModel, RampModel
from roman_datamodels.testing import utils as testutil
from roman_datamodels import datamodels

# Set parameters for multiple runs of data
args = "xstart, ystart, xsize, ysize, ngroups, instrument, exp_type"
test_data = [(1, 1, 2048, 2048, 2, 'WFI', 'WFI_IMAGE')]
ids = ["RampModel"]

@pytest.mark.parametrize(args, test_data, ids=ids)
def test_dq_im(xstart, ystart, xsize, ysize, ngroups, instrument, exp_type):
    """ Check that PIXELDQ is initialized with the information from the reference file.
    test that a flagged value in the reference file flags the PIXELDQ array"""

    csize = (ngroups, ysize, xsize)

    # create raw input data for step
    dm_ramp = testutil.mk_ramp(csize)
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    dq, dq_def = make_maskmodel_elements(ysize, xsize)

    # edit reference file with known bad pixel values
    dq[100, 100] = 2   # Dead pixel
    dq[200, 100] = 4   # Hot pixel
    dq[300, 100] = 8   # Unreliable_slope
    dq[400, 100] = 16  # NONLINEAR
    dq[500, 100] = 1   # Do_not_use
    dq[100, 200] = 3   # Dead pixel + do not use
    dq[200, 200] = 5   # Hot pixel + do not use
    dq[300, 200] = 9   # Unreliable slope + do not use
    dq[400, 200] = 17  # NONLINEAR + do not use

    # write mask model
    ref_data = testutil.mk_mask(arrays=csize[1:],dqsize=len(dq_def[0]))

    # Copy in maskmodel elemnts
    ref_data['dq'] = dq
    ref_data['dq_def'] = dq_def

    ref_data['meta']['instrument']['name'] = instrument

    # Instantiate model for pixel flag processing
    ref_data_model = MaskRefModel(ref_data)

    # run do_dqinit
    outfile = do_dqinit(dm_ramp, ref_data)

    dqdata = outfile.pixeldq

    # assert that the pixels read back in match the mapping from ref data to science data
    assert(dqdata[100, 100] == dqflags.pixel['DEAD'])
    assert(dqdata[200, 100] == dqflags.pixel['HOT'])
    assert(dqdata[300, 100] == dqflags.pixel['UNRELIABLE_SLOPE'])
    assert(dqdata[400, 100] == dqflags.pixel['NONLINEAR'])
    assert(dqdata[500, 100] == dqflags.pixel['DO_NOT_USE'])
    assert(dqdata[100, 200] == 1025)
    assert(dqdata[200, 200] == 2049)
    assert(dqdata[300, 200] == 16777217)
    assert (dqdata[400, 200] == 65537)


def test_groupdq():
    """Check that GROUPDQ extension is added to the data and all values are initialized to zero."""

    # size of integration
    instrument = 'WFI'
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = testutil.mk_ramp((ngroups, ysize, xsize))
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    dq, dq_def = make_maskmodel_elements(ysize, xsize)

    # write mask model
    ref_data = testutil.mk_mask(arrays=(ysize, xsize), dqsize=len(dq_def[0]))
    ref_data['meta']['instrument']['name'] = instrument

    # run the correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that GROUPDQ was created and initialized to zero
    groupdq = outfile.groupdq

    np.testing.assert_array_equal(np.full((ngroups, ysize, xsize), 0, dtype=int),
                                  groupdq, err_msg='groupdq not initialized to zero')


def test_err():
    """Check that a 4-D ERR array is initialized and all values are zero."""

    # size of integration
    instrument = 'WFI'
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = testutil.mk_ramp((ngroups, ysize, xsize))
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    dq, dq_def = make_maskmodel_elements(ysize, xsize)

    # write mask model
    ref_data = testutil.mk_mask(arrays=(ysize, xsize), dqsize=len(dq_def[0]))
    ref_data['meta']['instrument']['name'] = instrument

    # Filter out validation warnings from ref_data
    warnings.filterwarnings("ignore", category=ValidationWarning)

    # run correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that ERR array was created and initialized to zero
    errarr = outfile.err

    assert(errarr.ndim == 2)  # check that output err array is 2-D
    assert(np.all(errarr == 0))  # check that values are 0


def test_dq_add1_groupdq():
    """
    Test if the dq_init code set the groupdq flag on the first
    group to 'do_not_use' by adding 1 to the flag, not overwriting to 1
    Also test whether two flags on the same pixel are added together.
    """

    # size of integration
    instrument = 'WFI'
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = testutil.mk_ramp((ngroups, ysize, xsize))
    dm_ramp.meta.instrument.name = instrument

    # create a MaskModel elements for the dq input mask
    dq, dq_def = make_maskmodel_elements(ysize, xsize)

    # write reference file with known bad pixel values
    dq[505, 505] = 1   # Do_not_use
    dq[400, 500] = 3  # do_not_use and dead pixel

    # write mask model
    ref_data = testutil.mk_mask(arrays=(ysize, xsize), dqsize=len(dq_def[0]))
    ref_data['meta']['instrument']['name'] = instrument

    # Copy in maskmodel elemnts
    ref_data['dq'] = dq
    ref_data['dq_def'] = dq_def

    # Instantiate model for pixel flag processing
    ref_data_model = MaskRefModel(ref_data)

    # set a flag in the pixel dq
    dm_ramp.pixeldq[505, 505] = 4

    # run correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # test if pixels in pixeldq were incremented in value by 1
    assert(outfile.pixeldq[505, 505] == 5)  # check that previous dq flag is added to mask value
    assert(outfile.pixeldq[400, 500] == 1025)  # check two flags propagate correctly


# Commented out until CRDS updated with mask information
# # Set parameters for multiple runs of guider data
# # args = "xstart, ystart, xsize, ysize, nints, ngroups, instrument, exp_type, detector"
# # test_data = [(1, 1, 2048, 2048, 2, 2, 'FGS', 'FGS_ID-IMAGE', 'GUIDER1'),
# #              (1, 1, 1032, 1024, 1, 5, 'MIRI', 'MIR_IMAGE', 'MIRIMAGE')]
# # ids = ["GuiderRawModel-Image", "RampModel"]
# # Set parameters for multiple runs of data
# args = "xstart, ystart, xsize, ysize, ngroups, instrument, exp_type, detector"
# test_data = [(1, 1, 2048, 2048, 2, 'WFI', 'WFI_IMAGE','WFI04')]
# ids = ["RampModel"]
#
#
# @pytest.mark.parametrize(args, test_data, ids=ids)
# #def test_fullstep(xstart, ystart, xsize, ysize, nints, ngroups, instrument, exp_type, detector):
# def test_fullstep(xstart, ystart, xsize, ysize, ngroups, instrument, exp_type, detector):
#     """Test that the full step runs"""
#
#     # create raw input data for step
#     #dm_ramp = make_rawramp(instrument, ngroups, ysize, xsize, ystart, xstart, exp_type)
#     dm_ramp = testutil.mk_ramp((ngroups, ysize, xsize))
#
#     dm_ramp.meta.instrument.name = instrument
#     dm_ramp.meta.instrument.detector = detector
#     # dm_ramp.meta.observation.date = '2000-06-01'
#     # dm_ramp.meta.observation.time = '00:00:00'
#
#     # Instantiate model for pixel flag processing
#     dm_ramp_model = RampModel(dm_ramp)
#
#     # run the full step
#     #outfile = DQInitStep.call(dm_ramp)
#     outfile = DQInitStep.call(dm_ramp_model)
#
#     # test that a pixeldq frame has been initlialized
#     if instrument == "FGS":
#         assert outfile.dq.ndim == 2
#     else:
#         assert outfile.pixeldq.ndim == 2  # a 2-d pixeldq frame exists


def make_maskmodel_elements(ysize, xsize):
    # create a mask model for the dq_init step
    csize = (ysize, xsize)
    dq = np.zeros(csize, dtype=np.uint32)

    # define a dq_def extension
    dqdef = [(0, 1, 'DO_NOT_USE', 'Bad Pixel do not use'),
             (1, 2, 'DEAD', 'Dead Pixel'),
             (2, 4, 'HOT', 'Hot pixel'),
             (3, 8, 'UNRELIABLE_SLOPE', 'Large slope variance'),
             (4, 16, 'NONLINEAR', 'Pixel highly nonlinear'),
             (5, 32, 'REFERENCE_PIXEL', 'Reference Pixel')]

    dq_def = np.array((dqdef), dtype=table_definitions.DQ_DEF_DTYPE)

    return dq, dq_def
