import pytest
import roman_datamodels as rdm

from romancal.pipeline.mosaic_pipeline import MosaicPipeline

from .regtestdata import compare_asdf


@pytest.mark.bigdata
@pytest.mark.soctests
@pytest.fixture(scope="module")
def run_mos(rtdata_module):
    rtdata = rtdata_module

    # Test Pipeline
    rtdata.get_asn("WFI/image/L3_mosaic_asn.json")
    output = "r0099101001001001001_r274dp63x31y81_prompt_F158_i2d.asdf"
    rtdata.output = output
    args = [
        "roman_mos",
        rtdata.input,
    ]
    MosaicPipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_mos):
    return run_mos.output


@pytest.fixture(scope="module")
def output_model(output_filename):
    with rdm.open(output_filename) as model:
        yield model


@pytest.fixture(scope="module")
def truth_filename(run_mos):
    return run_mos.truth


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_resample_ran(output_model):
    # DMS373 test output is resampled to a skycell
    # FIXME this doesn't test if the output is a skyceell
    assert output_model.meta.cal_step.resample == "COMPLETE"
