"""
Purpose: record/query_step_status follow Roman cal_step_flag conventions and
tolerate missing meta.cal_step.
"""

import pytest
from roman_datamodels.datamodels import ImageModel

from romancal.datamodels.library import ModelLibrary
from romancal.flux.flux_step import apply_flux_correction
from romancal.stpipe import RomanStep
from romancal.stpipe.utilities import (
    ALLOWED_STEP_STATUSES,
    query_step_status,
    record_step_status,
)


def _image_without_cal_step():
    model = ImageModel.create_fake_data()
    if "cal_step" in model.meta:
        del model.meta["cal_step"]
    assert getattr(model.meta, "cal_step", None) is None
    return model


def test_query_missing_cal_step_returns_incomplete():
    """Purpose: missing meta.cal_step queries as INCOMPLETE, not an exception."""
    model = _image_without_cal_step()
    assert query_step_status(model, "flux") == "INCOMPLETE"


def test_record_creates_cal_step_and_sets_status():
    """Purpose: recording status creates meta.cal_step when absent."""
    model = _image_without_cal_step()
    record_step_status(model, "flux")
    assert "cal_step" in model.meta
    assert model.meta.cal_step.flux == "COMPLETE"
    assert query_step_status(model, "flux") == "COMPLETE"


@pytest.mark.parametrize("status", sorted(ALLOWED_STEP_STATUSES))
def test_record_accepts_all_roman_statuses(status):
    """Purpose: every rad cal_step_flag value can be recorded."""
    model = _image_without_cal_step()
    record_step_status(model, "flat_field", status=status)
    assert query_step_status(model, "flat_field") == status


def test_record_rejects_invalid_status():
    """Purpose: non-Roman status strings are rejected."""
    model = _image_without_cal_step()
    with pytest.raises(ValueError, match="Invalid cal_step status"):
        record_step_status(model, "flux", status="DONE")


def test_roman_step_record_step_status_complete_default():
    """Purpose: RomanStep.record_step_status defaults to COMPLETE on success."""
    model = _image_without_cal_step()
    step = RomanStep()
    step.record_step_status(model, "photom")
    assert query_step_status(model, "photom") == "COMPLETE"


def test_roman_step_record_step_status_requires_status_when_unsuccessful():
    """
    Purpose: unsuccessful steps must pass an explicit Roman status
    (SKIPPED/FAILED/N/A), not a JWST-style implicit FAILED.
    """
    model = _image_without_cal_step()
    step = RomanStep()
    with pytest.raises(ValueError, match="status must be set explicitly"):
        step.record_step_status(model, "tweakreg", success=False)
    step.record_step_status(model, "tweakreg", success=False, status="FAILED")
    assert query_step_status(model, "tweakreg") == "FAILED"
    step.record_step_status(model, "source_catalog", status="SKIPPED")
    assert query_step_status(model, "source_catalog") == "SKIPPED"


def test_flux_correction_tolerates_missing_cal_step():
    """
    Purpose: flux correction (issue #2213 example) no longer assumes cal_step
    exists; missing status is treated as not yet COMPLETE.
    """
    model = _image_without_cal_step()
    model.meta.photometry.conversion_megajanskys = 2.0
    original = model.data.copy()
    apply_flux_correction(model)
    assert model.data == pytest.approx(original * 2.0)


def test_record_and_query_model_library():
    """Purpose: ModelLibrary inputs are supported for record and query."""
    model = _image_without_cal_step()
    library = ModelLibrary([model])
    assert query_step_status(library, "flux") == "INCOMPLETE"
    record_step_status(library, "flux", status="COMPLETE")
    assert query_step_status(library, "flux") == "COMPLETE"
