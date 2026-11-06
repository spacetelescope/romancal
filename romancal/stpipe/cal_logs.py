"""Helpers for attaching calibration logs to Roman products.

Product policy (see https://github.com/spacetelescope/romancal/issues/2510):

* A product should carry logs relevant to producing that product.
* Full logs from earlier products should not be carried forward just because
  metadata is copied from an upstream model.
* L2 ELP association members should ideally hold only that member's logs.
* L3 mosaics keep input L2 logs on ``individual_image_cal_logs`` (not in
  ``meta.cal_logs``). Multiband catalogs must not embed L3 ``cal_logs`` in
  ``meta.image_metas`` (those logs already live on the L3 products).
"""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping, Sequence
from copy import deepcopy
from typing import Any

# Keys that must not be copied when seeding a new product's metadata from an
# upstream image/mosaic product.
_CAL_LOG_META_KEYS = frozenset({"cal_logs"})

# Keys omitted when building per-image metadata blobs for multiband catalogs.
# ``cal_logs`` is excluded so L3 production logs do not land in catalog parquet
# metadata via ``meta.image_metas``.
IMAGE_META_SKIP_KEYS = frozenset({"wcs", "individual_image_meta", "cal_logs"})


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
    Write calibration log records onto a product model.

    Prefers ``model.meta.cal_logs`` (current Image/Mosaic schema location).
    Falls back to top-level ``model.cal_logs`` only if meta cannot store the
    attribute.

    Parameters
    ----------
    model : DataModel
        Output product model.
    records : sequence of str or None
        Log lines to store. ``None`` is treated as an empty list.
    """
    payload = list(records) if records is not None else []
    meta = getattr(model, "meta", None)
    if meta is not None:
        if hasattr(meta, "__setitem__"):
            meta["cal_logs"] = payload
            return
        meta.cal_logs = payload
        return
    # Legacy / alternate placement used by some older mosaic products.
    if hasattr(model, "__setitem__"):
        model["cal_logs"] = payload
    else:
        model.cal_logs = payload


def strip_inherited_cal_logs(model_or_meta) -> None:
    """
    Remove upstream ``cal_logs`` from a model or meta mapping in place.

    Safe no-op when ``cal_logs`` is absent.
    """
    meta = _as_meta(model_or_meta)
    if meta is None:
        # Maybe a bare meta dict was not found; try model top-level.
        _pop_cal_logs(model_or_meta)
        return
    _pop_cal_logs(meta)
    # Also clear legacy top-level cal_logs on models.
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


def extract_image_meta(image_meta: Mapping[str, Any]) -> dict[str, Any]:
    """
    Copy per-image metadata for multiband ``image_metas``, excluding logs/WCS.

    Parameters
    ----------
    image_meta : mapping
        ``image_model.meta`` contents.

    Returns
    -------
    dict
        Deep-copied metadata suitable for ``meta.image_metas`` entries.
    """
    return {
        key: deepcopy(value)
        for key, value in image_meta.items()
        if key not in IMAGE_META_SKIP_KEYS
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
        # Prefer explicit meta attribute when present (DataModel).
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
