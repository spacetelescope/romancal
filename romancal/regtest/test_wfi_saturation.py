"""Tests for the saturation step"""

import os

import pytest
import roman_datamodels as rdm

from romancal.step import SaturationStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_saturation_image_step(rtdata, ignore_asdf_paths, resource_tracker, request):
    """Testing retrieval of best ref file for image data,
    and creation of a ramp file with CRDS selected saturation file applied."""

    input_file = "r0000101001001001001_0001_wfi01_f158_dqinit.asdf"
    rtdata.get_data(f"WFI/image/{input_file}")
    rtdata.input = input_file

    # Test CRDS
    step = SaturationStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "saturation")
    ref_file_name = os.path.split(ref_file_path)[-1]

    # Confirm that a WFI saturation reference file name was returned
    assert "roman_wfi_saturation" in ref_file_name

    # Test SaturationStep
    output = "r0000101001001001001_0001_wfi01_f158_saturation.asdf"
    rtdata.output = output

    args = ["romancal.step.SaturationStep", rtdata.input]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    ramp_out = rdm.open(rtdata.output)
    assert "roman.pixeldq" in ramp_out.to_flat_dict()

    rtdata.get_truth(f"truth/WFI/image/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_saturation_grism_step(rtdata, ignore_asdf_paths, resource_tracker, request):
    """Testing retrieval of best ref file for grism data,
    and creation of a ramp file with CRDS selected saturation file applied."""

    input_file = "r0000201001001001001_0001_wfi01_grism_dqinit.asdf"
    rtdata.get_data(f"WFI/grism/{input_file}")
    rtdata.input = input_file

    # Test CRDS
    step = SaturationStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "saturation")
    ref_file_name = os.path.split(ref_file_path)[-1]

    # Confirm that a WFI saturation reference file name was returned
    assert "roman_wfi_saturation" in ref_file_name

    # Test SaturationStep
    output = "r0000201001001001001_0001_wfi01_grism_saturation.asdf"
    rtdata.output = output

    args = ["romancal.step.SaturationStep", rtdata.input]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)
    resource_tracker.log(request)

    ramp_out = rdm.open(rtdata.output)
    assert "roman.pixeldq" in ramp_out.to_flat_dict()

    rtdata.get_truth(f"truth/WFI/grism/{output}")
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()
