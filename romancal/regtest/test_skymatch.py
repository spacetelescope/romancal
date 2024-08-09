import os
import pytest
from roman_datamodels import datamodels as rdm
from metrics_logger.decorators import metrics_logger
from romancal.datamodels.library import ModelLibrary
from romancal.pipeline.mosaic_pipeline import MosaicPipeline
from romancal.skymatch import SkyMatchStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


def is_not_none_or_empty(value):
    """
    Check if a value is neither None nor empty.

    Parameters
    ----------
    value : any
        The value to check. It can be of any type.

    Returns
    -------
    bool
        True if the value is neither None nor empty, False otherwise.

    Examples
    --------
    >>> is_not_none_or_empty(None)
    False
    >>> is_not_none_or_empty("")
    False
    >>> is_not_none_or_empty([])
    False
    >>> is_not_none_or_empty({})
    False
    >>> is_not_none_or_empty("Hello")
    True
    >>> is_not_none_or_empty([1, 2, 3])
    True
    >>> is_not_none_or_empty({"key": "value"})
    True
    """
    if value is None:
        return False
    if isinstance(value, (str, list, dict, set, tuple)) and len(value) == 0:
        return False
    return True


@metrics_logger()
@pytest.mark.bigdata
def test_skymatch_step(rtdata):
    """Test for the sky match step using imaging data."""

    rtdata.get_asn("WFI/image/L3_regtest_asn.json")
    rtdata.output = "skymatch_test"

    args = [
        "romancal.step.SkyMatchStep",
        rtdata.input,
        f"--output_file='{rtdata.output}'",
        "--subtract=True",
    ]
    RomanStep.from_cmdline(args)

    step = SkyMatchStep()
    suffix = step.default_suffix()

    n_members = len(rtdata.asn["products"][0]["members"])
    # create a list with the output filenames
    output_filenames = [f"{rtdata.output}_{i}_{suffix}.asdf" for i in range(n_members)]

    # instantiate a ModelLibrary object
    library = ModelLibrary(output_filenames)
    with library:
        for index, model in enumerate(library):
            # Perform DMS tests
            step.log.info(
                "DMSXXX MSG: SkyMatchStep added meta.background? :"
                f'  {hasattr(model.meta, "background")}'
            )
            assert hasattr(model.meta, "background")

            step.log.info(
                "DMSXXX MSG: SkyMatchStep populated meta.background? :"
                f"  {all(is_not_none_or_empty(v) for v in model.meta.background.values())}"
            )
            assert all(is_not_none_or_empty(v) for v in model.meta.background.values())

            library.shelve(model, index)


@metrics_logger()
@pytest.mark.bigdata
def test_l3_output_contains_background_info(rtdata):
    """Test for the presence of meta.background in L3 files."""

    rtdata.get_asn("WFI/image/L3_regtest_asn.json")
    rtdata.output = "test_background_info"

    args = [
        "roman_mos",
        rtdata.input,
        f"--output_file='{rtdata.output}'",
    ]
    MosaicPipeline.from_cmdline(args)

    step = MosaicPipeline()
    suffix = step.default_suffix()

    n_members = len(rtdata.asn["products"][0]["members"])
    # create a list with the output filenames
    output_filenames = [f"{rtdata.output}_{i}_{suffix}.asdf" for i in range(n_members)]

    # instantiate a ModelLibrary object
    library = ModelLibrary(output_filenames)
    with library:
        for index, model in enumerate(library):
            # Perform DMS tests
            step.log.info(
                "DMSXXX MSG: SkyMatchStep added meta.background? :"
                f'  {hasattr(model.meta, "background")}'
            )
            assert hasattr(model.meta, "background")

            step.log.info(
                "DMSXXX MSG: SkyMatchStep populated meta.background? :"
                f"  {all(is_not_none_or_empty(v) for v in model.meta.background.values())}"
            )
            assert all(is_not_none_or_empty(v) for v in model.meta.background.values())

            library.shelve(model, index)
