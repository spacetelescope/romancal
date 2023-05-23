from romancal.refpix.data import StandardView
from romancal.refpix.refpix import Control, run_steps

from . import reference_utils


def test_run_steps_regression(datamodel, ref_pix_ref):
    regression = StandardView.from_datamodel(datamodel).data.copy()
    regression_out = reference_utils.apply_correction(
        regression, ref_pix_ref.alpha, ref_pix_ref.gamma, ref_pix_ref.zeta
    )

    control = Control()
    result = run_steps(datamodel, ref_pix_ref, control)

    assert (result.data.value == regression_out).all()
