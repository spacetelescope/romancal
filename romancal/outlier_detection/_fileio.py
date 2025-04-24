import logging

import asdf

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_median(median_data, median_wcs, make_output_path):
    _save_intermediate_output(
        median_data,
        median_wcs,
        make_output_path("drizzled_median.asdf", suffix="median"),
    )


def save_drizzled(drizzled_model, make_output_path):
    input_path = drizzled_model.meta.filename.replace("_outlier_", "_")
    _save_intermediate_output(
        drizzled_model.data,
        drizzled_model.meta.wcs,
        make_output_path(input_path, suffix="outlier_coadd"),
    )


def _save_intermediate_output(data, wcs, output_path):
    asdf.AsdfFile({"data": data, "wcs": wcs}).write_to(output_path)
    log.info(f"Saved intermediate product to {output_path}")
