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

    Min/max coadd fields (``time_first``, ``time_last``,
    ``max_exposure_time``) and ordinary disagree-blend fields (including
    ``instrument.optical_element``) are updated here. Mean coadd fields are
    only accumulated; call :func:`finalize_catalog_metadata` after all inputs
    have been blended to set those finals.

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
            # max-like fields stay in the blender (same pattern as time_last)
            max_exptime = value.get("max_exposure_time")
            if max_exptime is not None:
                current = cat_model.meta[key].get("max_exposure_time")
                cat_model.meta[key]["max_exposure_time"] = (
                    float(max_exptime)
                    if current is None
                    else float(max(current, max_exptime))
                )
            # means are awkward as running blends; accumulate for finalize
            if value.get("time_mean") is not None:
                time_means.append(value["time_mean"])
            if value.get("exposure_time") is not None:
                exposure_times.append(value["exposure_time"])
        else:
            # set non-matching metadata values to None (covers optical_element)
            for subkey, subvalue in value.items():
                if cat_model.meta[key].get(subkey, None) != subvalue:
                    cat_model.meta[key][subkey] = None


def finalize_catalog_metadata(cat_model, time_means, exposure_times):
    """
    Finalize blended catalog metadata after all L3 inputs are processed.

    Sets coadd_info mean fields from values accumulated during blending.
    Min/max coadd fields and multi-filter optical_element nulling are handled
    in :func:`blend_image_metadata`.

    Parameters
    ----------
    cat_model : MultibandSourceCatalogModel
        The multiband catalog model being built (modified in place).

    time_means : list
        Accumulated mean observation times from input L3 coadds.

    exposure_times : list
        Accumulated exposure times from input L3 coadds.
    """
    if time_means:
        cat_model.meta.coadd_info.time_mean = Time(time_means).mean()
    if exposure_times:
        cat_model.meta.coadd_info.exposure_time = float(np.mean(exposure_times))
