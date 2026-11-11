"""Helpers for reading/writing product calibration log lists.

Used mainly by ExposurePipeline when slicing multi-member ELP log streams,
and by catalog construction when seeding metadata without inheriting
upstream ``cal_logs``.
"""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping, Sequence
from copy import deepcopy
from typing import Any

# Keys that must not be copied when seeding a new product's metadata from an
# upstream image/mosaic product.
_CAL_LOG_META_KEYS = frozenset({"cal_logs"})


def get_cal_logs(model_or_meta) -> list[str]:
    """
    Return cal_logs from a model or meta mapping as a plain list.

    Parameters
    ----------
    model_or_meta : DataModel or mapping
        Object with ``meta.cal_logs`` or a meta mapping that may contain
        ``cal_logs``.

    Returns
    -------
    list of str
        Copy of the log records, or an empty list if none are present.
    """
    meta = _as_meta(model_or_meta)
    if meta is None:
        return []
    logs = meta.get("cal_logs") if hasattr(meta, "get") else None
    if logs is None and hasattr(meta, "cal_logs"):
        logs = meta.cal_logs
    if not logs:
        return []
    return list(logs)


def set_product_cal_logs(model, records: Sequence[str] | None) -> None:
    """
    Write calibration log records onto ``model.meta.cal_logs``.

    This is equivalent to ``model.meta.cal_logs = list(records)`` for normal
    Image/Mosaic models. The helper exists so ELP member-scoping and tests can
    assign log lists without repeating that assignment.

    Parameters
    ----------
    model : DataModel
        Output product model.
    records : sequence of str or None
        Log lines to store. ``None`` is treated as an empty list.
    """
    payload = list(records) if records is not None else []
    model.meta["cal_logs"] = payload


def strip_inherited_cal_logs(model_or_meta) -> None:
    """
    Remove upstream ``cal_logs`` from a model or meta mapping in place.

    Safe no-op when ``cal_logs`` is absent.
    """
    meta = _as_meta(model_or_meta)
    if meta is None:
        _pop_cal_logs(model_or_meta)
        return
    _pop_cal_logs(meta)
    if model_or_meta is not meta:
        _pop_cal_logs(model_or_meta)


def copy_meta_without_cal_logs(input_meta) -> dict[str, Any]:
    """
    Deep-copy product metadata while dropping calibration logs.

    Parameters
    ----------
    input_meta : mapping
        Source ``model.meta`` (or equivalent).

    Returns
    -------
    dict
        Independent metadata copy without ``cal_logs``.
    """
    if input_meta is None:
        return {}
    if hasattr(input_meta, "items"):
        items = input_meta.items()
    else:
        items = dict(input_meta).items()
    return {
        key: deepcopy(value) for key, value in items if key not in _CAL_LOG_META_KEYS
    }


def slice_log_records(
    records: Sequence[str], start: int, end: int | None = None
) -> list[str]:
    """
    Return ``records[start:end]`` as a new list.

    Parameters
    ----------
    records : sequence of str
        Accumulated step/pipeline log records.
    start : int
        Inclusive start index.
    end : int or None
        Exclusive end index; defaults to ``len(records)``.
    """
    if end is None:
        end = len(records)
    start = max(0, start)
    end = max(start, end)
    return list(records[start:end])


def _as_meta(model_or_meta):
    if model_or_meta is None:
        return None
    if isinstance(model_or_meta, Mapping) or hasattr(model_or_meta, "get"):
        meta = getattr(model_or_meta, "meta", None)
        if meta is not None and meta is not model_or_meta:
            return meta
        return model_or_meta
    return getattr(model_or_meta, "meta", None)


def _pop_cal_logs(obj) -> None:
    if obj is None:
        return
    if isinstance(obj, MutableMapping) or hasattr(obj, "pop"):
        if "cal_logs" in obj:
            obj.pop("cal_logs", None)
        return
    if hasattr(obj, "__contains__") and hasattr(obj, "__delitem__"):
        if "cal_logs" in obj:
            del obj["cal_logs"]
        return
    if hasattr(obj, "cal_logs"):
        obj.cal_logs = []
