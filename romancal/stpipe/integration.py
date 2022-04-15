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
        ("romancal.pipeline.ExposurePipeline", 'roman_elp', True),
        ("romancal.step.DarkCurrentStep", None, False),
        ("romancal.step.DQInitStep", None, False),
        ("romancal.step.FlatFieldStep", None, False),
        ("romancal.step.JumpStep", None, False),
        ("romancal.step.LinearityStep", None, False),
        ("romancal.step.PhotomStep", None, False),
        ("romancal.step.RampFitStep", None, False),
        ("romancal.step.SaturationStep", None, False),
        ("romancal.step.AssignWcsStep", None, False),
    ]
