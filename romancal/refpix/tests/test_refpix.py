import numpy as np
from numpy.testing import assert_allclose

from romancal.refpix._data import StandardView
from romancal.refpix.refpix import run_steps

from . import reference_utils


def test_run_steps_regression(datamodel, ref_pix_ref):
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
