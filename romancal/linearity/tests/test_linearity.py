"""

Unit tests for linearity correction

"""

import numpy as np
import pytest
from astropy import units as u
from roman_datamodels import maker_utils
from roman_datamodels.datamodels import LinearityRefModel, ScienceRawModel

from romancal.dq_init import DQInitStep
from romancal.lib import dqflags
from romancal.linearity import LinearityStep


@pytest.mark.parametrize(
    "instrument, exptype",
    [
        ("WFI", "WFI_IMAGE"),
    ],
)
def test_linearity_coeff(instrument, exptype):
    """Test that the basic inferface works for data requiring a DQ reffile"""

    # Set test size
    shape = (5, 20, 20)

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

    result = DQInitStep.call(wfi_sci_raw_model)
    result = LinearityStep.call(result, override_linearity="N/A")

    assert result.meta.cal_step.linearity == "SKIPPED"

    # Set coefficient values in reference file to check the algorithm
    # Equation is DNcorr = L0 + L1*DN(i) + L2*DN(i)^2 + L3*DN(i)^3 + L4*DN(i)^4
    # DN(i) = signal in pixel, Ln = coefficient from ref file
    # L0 = 0 for all pixels for CDP6

    coeffs = np.zeros(shape, dtype=np.float32)
    coeffs[1, :, :].fill(0.85)
    coeffs[2, :, :].fill(4.62e-06)
    coeffs[3, :, :].fill(-6.16e-11)
    coeffs[4, :, :].fill(7.23e-16)
    # coeffs = np.tile(coeffs,(5,20,4))
    # set one set of coeffs to 0 to check that the data values are not altered.
    # and that the DQ flag is set to NO_LIN_CORR
    coeffs[0:, 5, 5] = 0.0
    # set one of the coeffs to NaN and check that the DQ flag is set to NO_LIN_CORR
    coeffs[0:, 6, 5] = np.nan
    # save the pixel values to make sure they are not altered
    pixels_55 = result.data[0:, 5, 5]
    pixels_65 = result.data[0:, 6, 5]
    linref = maker_utils.mk_linearity(shape=shape, coeffs=coeffs)
    linref_model = LinearityRefModel(linref)

    LinearityStep.call(result, override_linearity=linref_model)

    assert result.meta.cal_step.linearity == "COMPLETE"
    assert result.pixeldq[5, 5] == dqflags.pixel["NO_LIN_CORR"]
    assert result.pixeldq[6, 5] == dqflags.pixel["NO_LIN_CORR"]
    np.testing.assert_array_equal(result.data[0:, 5, 5], pixels_55)
    np.testing.assert_array_equal(result.data[0:, 6, 5], pixels_65)
