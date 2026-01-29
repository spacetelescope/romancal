"""

Unit tests for linearity correction

"""

import numpy as np
import pytest
from astropy.time import Time
from roman_datamodels import dqflags
from roman_datamodels.datamodels import (
    IntegralnonlinearityRefModel,
    InverselinearityRefModel,
    LinearityRefModel,
    ScienceRawModel,
)

from romancal.dq_init import DQInitStep
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
    wfi_sci_raw_model = ScienceRawModel.create_fake_data(shape=shape)
    wfi_sci_raw_model.meta.exposure.start_time = Time(
        "2024-01-03T00:00:00.0", format="isot", scale="utc"
    )
    wfi_sci_raw_model.meta.instrument.name = instrument
    wfi_sci_raw_model.meta.instrument.detector = "WFI01"
    wfi_sci_raw_model.meta.instrument.optical_element = "F158"
    wfi_sci_raw_model.meta["guide_star"]["window_xstart"] = 1012
    wfi_sci_raw_model.meta["guide_star"]["window_xsize"] = 16
    wfi_sci_raw_model.meta.exposure.type = exptype
    wfi_sci_raw_model.data = np.ones(shape, dtype=np.uint16)

    result = DQInitStep.call(wfi_sci_raw_model)
    result = LinearityStep.call(
        result,
        override_linearity="N/A",
        override_inverselinearity="N/A",
        override_integralnonlinearity="N/A",
    )

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
    linref_model = LinearityRefModel.create_fake_data(
        {"coeffs": coeffs}, shape=shape[1:]
    )
    ilinref_model = InverselinearityRefModel.create_fake_data(
        {"coeffs": coeffs.copy()}, shape=shape[1:]
    )

    result.meta.exposure.read_pattern = [[1], [2], [3], [4], [5]]

    LinearityStep.call(
        result,
        override_linearity=linref_model,
        override_inverselinearity=ilinref_model,
        override_integralnonlinearity="N/A",
    )

    assert result.meta.cal_step.linearity == "COMPLETE"
    assert result.pixeldq[5, 5] == dqflags.pixel["NO_LIN_CORR"]
    assert result.pixeldq[6, 5] == dqflags.pixel["NO_LIN_CORR"]
    np.testing.assert_array_equal(result.data[0:, 5, 5], pixels_55)
    np.testing.assert_array_equal(result.data[0:, 6, 5], pixels_65)


def test_linearity_correction_values():
    """Test linearity correction with simple coefficients and INL."""

    shape = (5, 20, 256)

    wfi_sci_raw_model = ScienceRawModel.create_fake_data(shape=shape)
    wfi_sci_raw_model.meta.exposure.start_time = Time(
        "2024-01-03T00:00:00.0", format="isot", scale="utc"
    )
    wfi_sci_raw_model.meta.instrument.name = "WFI"
    wfi_sci_raw_model.meta.instrument.detector = "WFI01"
    wfi_sci_raw_model.meta.instrument.optical_element = "F158"
    wfi_sci_raw_model.meta["guide_star"]["window_xstart"] = 1012
    wfi_sci_raw_model.meta["guide_star"]["window_xsize"] = 16
    wfi_sci_raw_model.meta.exposure.type = "WFI_IMAGE"
    wfi_sci_raw_model.data = np.zeros(shape, dtype=np.uint16)
    wfi_sci_raw_model.data[:, 0, :] = 0
    wfi_sci_raw_model.data[:, 1, :] = 100

    result = DQInitStep.call(wfi_sci_raw_model)

    # Simple linearity: output = 2 * input
    lin_coeffs = np.zeros(shape, dtype=np.float32)
    lin_coeffs[1, :, :] = 2.0

    # Simple inverse linearity: output = 0.5 * input
    ilin_coeffs = np.zeros(shape, dtype=np.float32)
    ilin_coeffs[1, :, :] = 0.5

    linref_model = LinearityRefModel.create_fake_data(
        {"coeffs": lin_coeffs}, shape=shape[1:]
    )
    ilinref_model = InverselinearityRefModel.create_fake_data(
        {"coeffs": ilin_coeffs}, shape=shape[1:]
    )

    # Create INL model with correction = 1 everywhere
    inl_table_data = {}
    for start_col in range(0, shape[2], 128):
        channel_num = start_col // 128 + 1
        attr_name = f"science_channel_{channel_num:02d}"
        inl_table_data[attr_name] = {"correction": np.ones(65536, dtype="f4")}

    inlref_model = IntegralnonlinearityRefModel.create_fake_data({
        "value": np.arange(65536, dtype="f4"),
        "inl_table": inl_table_data,
    })

    result.meta.exposure.read_pattern = [[1], [2], [3], [4], [5]]

    LinearityStep.call(
        result,
        override_linearity=linref_model,
        override_inverselinearity=ilinref_model,
        override_integralnonlinearity=inlref_model,
    )

    assert result.meta.cal_step.linearity == "COMPLETE"

    # the integral non-linearity is applied to some simulated
    # individual non-linear reads.  Then the classical linearity
    # correction is applied on top of that.
    # In this case the ramp slope is zero so the resultant ->
    # read simulation is trivial doesn't do anything,
    # so this just applies the INL and then the CNL.
    
    # integral non-linearity corrects 0 -> 1, classical non-linearity
    # doubles 1->2.
    np.testing.assert_array_equal(result.data[:, 0, :], 2)

    # integral non-linearity corrects 100 -> 101, classical non-linearity
    # doubles to 202.
    np.testing.assert_array_equal(result.data[:, 1, :], 202)
