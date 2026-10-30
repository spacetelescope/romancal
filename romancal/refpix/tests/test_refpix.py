import numpy as np
import pytest
from numpy.testing import assert_allclose

from romancal.refpix._data import StandardView
from romancal.refpix.refpix import run_steps

from . import reference_utils


@pytest.mark.parametrize("force_state", ["default", "left", "right", "none"])
def test_dynamic_weight_scaffolding(datamodel, ref_pix_ref, force_state):
    """
    Validates that the pipeline can accept and process the 4 dynamic
    weight states without broadcasting or dimensional errors.
    """

    # We need to ensure that there is a read pattern compatible with
    # the size of the data array.
    nreads = datamodel.data.shape[0]
    datamodel.meta.exposure.read_pattern = [[i] for i in range(1, nreads + 1)]
    datamodel.meta.exposure.frame_time = 1.0

    # Create a clean copy so we don't poison the fixture for subsequent tests
    test_model = datamodel.copy()

    processed_model = run_steps(test_model, ref_pix_ref, weight_mode=force_state)

    # Assert the pipeline successfully returned a RampModel
    assert processed_model is not None
    assert processed_model.data.shape == datamodel.data.shape


def test_run_steps_regression(datamodel, ref_pix_ref):
    # We need to ensure that there is a read pattern compatible with
    # the size of the data array.
    nreads = datamodel.data.shape[0]
    datamodel.meta.exposure.read_pattern = [[i] for i in range(1, nreads + 1)]
    datamodel.meta.exposure.frame_time = 1.0

    regression = StandardView.from_datamodel(datamodel).data.copy()
    regression_out = reference_utils.apply_correction(
        regression, ref_pix_ref.alpha, ref_pix_ref.gamma, ref_pix_ref.zeta
    )

    result = run_steps(
        datamodel,
        ref_pix_ref,
        remove_offset=True,
        remove_trends=True,
        cosine_interpolate=True,
        fft_interpolate=True,
    )

    # tolerance corresponds to one quantum of float32 precision for the largest values
    # and addresses different numerical approaches in different scipy versions
    assert_allclose(
        result.data,
        regression_out,
        rtol=1e-5,
        atol=1e-6 * np.abs(regression_out).max(),
    )
    # regression_out does not return amp33 data
