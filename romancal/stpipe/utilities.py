"""Utilities for Roman calibration step status metadata."""

from __future__ import annotations

from collections.abc import Sequence

from roman_datamodels.datamodels import DataModel

from romancal.datamodels.library import ModelLibrary

__all__ = ["ALLOWED_STEP_STATUSES", "query_step_status", "record_step_status"]

# Values allowed by rad enums/cal_step_flag (roman_datamodels schemas).
ALLOWED_STEP_STATUSES = frozenset(
    {
        "INCOMPLETE",
        "COMPLETE",
        "SKIPPED",
        "FAILED",
        "N/A",
    }
)

# Default when cal_step is missing or the step has never been recorded.
_NOT_SET = "INCOMPLETE"


def _ensure_cal_step(model: DataModel) -> None:
    """Create ``model.meta.cal_step`` if it does not exist."""
    if getattr(model.meta, "cal_step", None) is None:
        model.meta["cal_step"] = {}


def _as_model_list(models):
    """Normalize input into an iterable of datamodels without mutating libraries."""
    if isinstance(models, ModelLibrary):
        return models
    if isinstance(models, DataModel):
        return [models]
    if isinstance(models, Sequence) and not isinstance(models, str | bytes):
        return list(models)
    raise TypeError(
        "models must be a DataModel, sequence of DataModels, or ModelLibrary; "
        f"got {type(models)!r}"
    )


def _validate_status(status: str) -> str:
    if status not in ALLOWED_STEP_STATUSES:
        allowed = ", ".join(sorted(ALLOWED_STEP_STATUSES))
        raise ValueError(
            f"Invalid cal_step status {status!r}; expected one of: {allowed}"
        )
    return status


def record_step_status(models, cal_step, status="COMPLETE"):
    """
    Record completion status for a calibration step in model metadata.

    Ensures ``meta.cal_step`` exists before writing, so models that omit this
    optional attribute (for example older simulator products) do not raise.

    Status values follow the Roman ``cal_step_flag`` enum:
    ``INCOMPLETE``, ``COMPLETE``, ``SKIPPED``, ``FAILED``, and ``N/A``.

    Parameters
    ----------
    models : `~roman_datamodels.datamodels.DataModel`, sequence, or `~romancal.datamodels.ModelLibrary`
        One or more models to update.
    cal_step : str
        Name of the calibration step (for example ``"flux"`` or ``"flat_field"``).
    status : str, optional
        Step status string. Default is ``"COMPLETE"``.
    """
    status = _validate_status(status)
    models = _as_model_list(models)

    if isinstance(models, ModelLibrary):
        with models:
            for model in models:
                _ensure_cal_step(model)
                model.meta.cal_step[cal_step] = status
                models.shelve(model)
        return

    for model in models:
        _ensure_cal_step(model)
        model.meta.cal_step[cal_step] = status


def query_step_status(models, cal_step):
    """
    Query completion status for a calibration step from model metadata.

    Parameters
    ----------
    models : `~roman_datamodels.datamodels.DataModel`, sequence, or `~romancal.datamodels.ModelLibrary`
        Model or collection to query. For collections, the first model is used.
    cal_step : str
        Name of the calibration step.

    Returns
    -------
    str
        Recorded status, or ``"INCOMPLETE"`` if ``meta.cal_step`` or the step
        entry is missing (Roman default for an unset step).
    """
    models = _as_model_list(models)

    if isinstance(models, ModelLibrary):
        with models:
            model = models.borrow(0)
            models.shelve(model, modify=False)
    else:
        model = models[0]

    cal_step_meta = getattr(model.meta, "cal_step", None)
    if cal_step_meta is None:
        return _NOT_SET

    # Support both attribute and mapping access on DNode / dict-like objects.
    try:
        value = cal_step_meta[cal_step]
    except Exception:
        value = getattr(cal_step_meta, cal_step, _NOT_SET)

    if value is None:
        return _NOT_SET
    return value
