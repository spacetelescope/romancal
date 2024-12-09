"""
Entry point implementations.
"""


def get_steps():
    """
    Return tuples describing the stpipe.Step subclasses provided
    by this package.  This method is registered with the stpipe.steps
    entry point.

    Returns
    -------
    list of tuple (str, str, bool)
        The first element each tuple is a fully-qualified Step
        subclass name.  The second element is an optional class
        alias.  The third element indicates that the class
        is a subclass of Pipeline.
    """
    # Unit tests ensure that this list is kept in sync with the actual
    # class definitions.  We need to avoid importing romancal.pipeline and
    # romancal.step to keep the CLI snappy.
    return [
        ("romancal.pipeline.ExposurePipeline", "roman_elp", True),
        ("romancal.pipeline.MosaicPipeline", "roman_mos", True),
        ("romancal.step.DarkCurrentStep", "dark", False),
        ("romancal.step.DQInitStep", "dq_init", False),
        ("romancal.step.FlatFieldStep", "flat_field", False),
        ("romancal.step.FluxStep", "flux", False),
        ("romancal.step.LinearityStep", "linearity", False),
        ("romancal.step.PhotomStep", "photom", False),
        ("romancal.step.RampFitStep", "ramp_fit", False),
        ("romancal.step.RefPixStep", "refpix", False),
        ("romancal.step.SaturationStep", "saturation", False),
        ("romancal.step.AssignWcsStep", "assign_wcs", False),
        ("romancal.step.OutlierDetectionStep", "outlier_detection", False),
        ("romancal.step.SkyMatchStep", "skymatch", False),
        ("romancal.step.TweakRegStep", "tweakreg", False),
        ("romancal.step.ResampleStep", "resample", False),
        ("romancal.step.SourceCatalogStep", "source_catalog", False),
        ("romancal.step.MultibandCatalogStep", "multiband_catalog", False),
    ]
