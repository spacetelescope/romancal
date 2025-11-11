from copy import deepcopy

# Metadata keys to skip when accumulating image metadata
_SKIP_IMAGE_META_KEYS = {"wcs", "individual_image_meta"}

# Metadata keys to skip when blending metadata
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
    2. Appends it to the catalog's image_meta list
    3. Blends metadata values across filters, setting mismatches to None
    4. Handles special cases like coadd_info timing information
    5. Updates file_date to the earliest date

    This function modifies cat_model, time_means, and exposure_times in
    place.

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
            time_means.append(value["time_mean"])
            exposure_times.append(value["exposure_time"])
        else:
            # set non-matching metadata values to None
            for subkey, subvalue in value.items():
                if cat_model.meta[key].get(subkey, None) != subvalue:
                    cat_model.meta[key][subkey] = None
