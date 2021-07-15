import os
import pytest

from romancal.stpipe import RomanStep
from romancal.step import FlatFieldStep
import roman_datamodels as rdm
from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_flat_field_image_step(rtdata, ignore_asdf_paths):

    rtdata.get_data("WFI/image/l2_0001_rate.asdf")
    rtdata.input = "l2_0001_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "flat")
    ref_file_name = os.path.split(ref_file_path)[-1]
    assert "roman_wfi_flat" in ref_file_name

    # Test FlatFieldStep
    output = "l2_0001_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)

@pytest.mark.bigdata
def test_flat_field_grism_step(rtdata, ignore_asdf_paths):

    rtdata.get_data("WFI/grism/l2_0001_grism_rate.asdf")
    rtdata.input = "l2_0001_grism_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "flat")
    ref_file_name = os.path.split(ref_file_path)[-1]
    assert "roman_wfi_flat" in ref_file_name

    # Test FlatFieldStep
    output = "l2_0001_grism_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/grism/{output}")
    compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
