"""L3→L3 metadata blending for multiband source catalogs."""

from __future__ import annotations

from copy import deepcopy

import numpy as np
from astropy.time import Time

# Metadata keys to skip when accumulating image metadata
_SKIP_IMAGE_META_KEYS = {"wcs", "individual_image_meta"}

# Metadata keys to skip when blending metadata (keep first input)
_SKIP_BLEND_KEYS = {"wcsinfo"}


def blend_image_metadata(
    image_model,
    cat_model,
    time_means,
    exposure_times,
    max_exposure_times,
):
    """
    Accumulate and blend metadata from an individual filter image into
    the catalog.

    This function:
    1. Extracts relevant metadata from the input image model
    2. Appends it to the catalog's image_metas list
    3. Blends metadata values across filters, setting mismatches to None
    4. Handles special cases like coadd_info timing information
    5. Updates file_date to the earliest date

    This function modifies cat_model, time_means, exposure_times, and
    max_exposure_times in place. Call :func:`finalize_catalog_metadata`
    after all inputs have been blended to set mean/max coadd fields and
    multi-filter optical_element.

    Parameters
    ----------
    image_model : ImageModel or MosaicModel
        The input image model for a single filter.

    cat_model : MultibandSourceCatalogModel
        The multiband catalog model being built.

    time_means : list
        List to accumulate mean observation times (modified in place).

    exposure_times : list
        List to accumulate exposure times (modified in place).

    max_exposure_times : list
        List to accumulate max exposure times (modified in place).
    """
    # Accumulate image metadata
    image_meta = {
        k: deepcopy(v)
        for k, v in image_model["meta"].items()
        if k not in _SKIP_IMAGE_META_KEYS
    }
    cat_model.meta.image_metas.append(image_meta)

    # Blend model with catalog metadata
    if image_model.meta.file_date < cat_model.meta.image.file_date:
        cat_model.meta.image.file_date = image_model.meta.file_date

    for key, value in image_meta.items():
        if key in _SKIP_BLEND_KEYS:
            continue
        if not isinstance(value, dict):
            # skip blending of single top-level values
            continue
        if key not in cat_model.meta:
            # skip blending if the key is not in the catalog meta
            continue
        if key == "coadd_info":
            cat_model.meta[key]["time_first"] = min(
                cat_model.meta[key]["time_first"], value["time_first"]
            )
            cat_model.meta[key]["time_last"] = max(
                cat_model.meta[key]["time_last"], value["time_last"]
            )
            if value.get("time_mean") is not None:
                time_means.append(value["time_mean"])
            if value.get("exposure_time") is not None:
                exposure_times.append(value["exposure_time"])
            max_exptime = value.get("max_exposure_time")
            if max_exptime is not None:
                max_exposure_times.append(max_exptime)
        else:
            # set non-matching metadata values to None
            for subkey, subvalue in value.items():
                if cat_model.meta[key].get(subkey, None) != subvalue:
                    cat_model.meta[key][subkey] = None


def finalize_catalog_metadata(
    cat_model,
    time_means,
    exposure_times,
    max_exposure_times,
):
    """
    Finalize blended catalog metadata after all L3 inputs are processed.

    Sets coadd_info mean/max fields from accumulated values and forces
    top-level optical_element to None when multiple filters contribute.
    Per-filter optical elements remain available in image_metas.

    Parameters
    ----------
    cat_model : MultibandSourceCatalogModel
        The multiband catalog model being built (modified in place).

    time_means : list
        Accumulated mean observation times from input L3 coadds.

    exposure_times : list
        Accumulated exposure times from input L3 coadds.

    max_exposure_times : list
        Accumulated max exposure times from input L3 coadds.
    """
    if time_means:
        cat_model.meta.coadd_info.time_mean = Time(time_means).mean()
    if exposure_times:
        cat_model.meta.coadd_info.exposure_time = float(np.mean(exposure_times))
    if max_exposure_times:
        cat_model.meta.coadd_info.max_exposure_time = float(np.max(max_exposure_times))

    # Explicit multi-filter rule: top-level optical_element cannot hold a
    # list (enum), so set None when more than one distinct filter is present.
    # Individual filters remain in image_metas.
    filters = []
    for image_meta in cat_model.meta.get("image_metas", []):
        instrument = image_meta.get("instrument") or {}
        optical_element = instrument.get("optical_element")
        if optical_element is not None:
            filters.append(optical_element)

    if len(set(filters)) > 1:
        cat_model.meta.instrument.optical_element = None
