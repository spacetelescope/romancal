import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_median(example_model, median_data, median_wcs, make_output_path):
    _save_intermediate_output(
        _make_median_model(example_model, median_data, median_wcs),
        "median",
        make_output_path,
    )


def save_drizzled(drizzled_model, make_output_path):
    _save_intermediate_output(drizzled_model, "outlier_coadd", make_output_path)


def _make_median_model(example_model, data, wcs):
    model = example_model.copy()
    model.data = data
    model.meta.filename = "drizzled_median.asdf"
    model.meta.wcs = wcs
    return model


def _save_intermediate_output(model, suffix, make_output_path):
    """
    Ensure all intermediate outputs from OutlierDetectionStep have consistent file naming conventions

    Notes
    -----
    self.make_output_path() is updated globally for the step in the main pipeline
    to include the asn_id in the output path, so no need to handle it here.
    """

    # outlier_?2d is not a known suffix, and make_output_path cannot handle an
    # underscore in an unknown suffix, so do a manual string replacement
    input_path = model.meta.filename.replace("_outlier_", "_")

    output_path = make_output_path(input_path, suffix=suffix)
    model.save(output_path)
    log.info(f"Saved {suffix} model in {output_path}")
