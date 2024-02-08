from astropy import units as u

from romancal.refpix.data import StandardView
from romancal.refpix.refpix import run_steps

from . import reference_utils


def test_run_steps_regression(datamodel, ref_pix_ref):
    regression = StandardView.from_datamodel(datamodel).data.copy()
    regression_out = reference_utils.apply_correction(
        regression, ref_pix_ref.alpha, ref_pix_ref.gamma, ref_pix_ref.zeta
    )

    result = run_steps(datamodel, ref_pix_ref, True, True, True, True)

    assert (result.data.value == regression_out).all()
    # regression_out does not return amp33 data

    # Check the units
    assert result.data.unit == u.DN
    assert result.amp33.unit == u.DN
