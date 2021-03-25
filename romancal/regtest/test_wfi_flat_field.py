import os
import pytest
from asdf.commands import diff as asdf_diff

from romancal.stpipe import RomanStep
from romancal import datamodels


@pytest.fixture(scope='module')
def run_flat_field(rtdata_module):
    """ Run flat field on Level 2 imaging data."""

    rtdata = rtdata_module
    rtdata.get_data("WFI/image/l2_0001_rate.asdf")
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    RomanStep.from_cmdline(args)


@pytest.mark.bigdata
def test_flat_field_step(run_flat_field, rtdata_module):
    rtdata = rtdata_module
    rtdata.input = "l2_0001_rate.asdf"
    output = "l2_0001_rate_flatfieldstep.asdf"
    rtdata.output = output

    rtdata.get_truth(f"truth/WFI/image/{output}")

    report = asdf_diff([rtdata.output, rtdata.truth], False)
    assert report is not None, report


@pytest.fixture(scope='module')
def test_crds_retrieval(_jail, rtdata_module):
    """ Test retrieving a flat file from CRDS."""

    rtdata = rtdata_module
    rtdata.get_data("WFI/image/level2.asdf")
    step = FlatFieldStep()
    model = datamodels.ImageModel(rtdata.input)
    ref_file_path = step.get_reference_file(model, "flat")
    ref_file_name = os.path.split(ref_file_path)[-1]
    assert ref_file_name == "roman_wfi_flat_0003.asdf"
