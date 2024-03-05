""" Regression tests for the High Level Processing Pipeline"""
import pytest

from romancal.pipeline.highlevel_pipeline import HighLevelPipeline
from romancal.stpipe.core import RomanStep


@pytest.mark.xfail(reason='See RCAL-777')
@pytest.mark.bigdata
def test_flux(hlp):
    pass


# ########
# Fixtures
# ########
@pytest.fixture(scope='module')
def hlp(rtdata_module):
    """Execute the HighLevelPipeline"""
    rtdata = rtdata_module

    # Setup input data for an association
    input_data = [
        "r0000101001001001001_01101_0001_WFI01_cal.asdf",
        "r0000101001001001001_01101_0002_WFI01_cal.asdf",
    ]
    [rtdata.get_data(f"WFI/image/{data}") for data in input_data]
    asnfn = "mosaic_asn.json"
    rtdata.get_data(f"WFI/image/{asnfn}")

    # Execute pipeline
    rtdata.output = 'result_hlp.asdf'

    args = [
        'roman_hlp',
        rtdata.input,
        f'--output_file={rtdata.output}',
        '--steps.outlier_detection.skip=true',
    ]

    RomanStep.from_cmdline(args)

    return rtdata.output
