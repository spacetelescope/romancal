import os
import pytest

from romancal.stpipe import RomanStep
from romancal.step import SaturationStep
import roman_datamodels as rdm
from .regtestdata import compare_asdf

@pytest.mark.bigdata
def test_saturation_image_step(rtdata, ignore_asdf_paths):
    # Testing retrieval of best ref file for image data,
    # and creation of a ramp file with CRDS selected saturation file applied.

    rtdata.get_data("WFI/image/ramp.asdf")
    rtdata.input = "ramp.asdf"

    # Test CRDS
    step = SaturationStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "saturation")
    ref_file_name = os.path.split(ref_file_path)[-1]

    # Confirm that a WFI saturation reference file name was returned
    assert "roman_wfi_saturation" in ref_file_name

    # Test SaturationStep
    output = "ramp_saturationstep.asdf"
    rtdata.output = output

    args = ["romancal.step.SaturationStep", rtdata.input]
    RomanStep.from_cmdline(args)

    ramp_out = rdm.open(rtdata.output)
    assert ("roman.pixeldq" in ramp_out.to_flat_dict())

    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)


@pytest.mark.bigdata
def test_saturation_grism_step(rtdata, ignore_asdf_paths):
    # Testing retrieval of best ref file for grism data,
    # and creation of a ramp file with CRDS selected saturation file applied.

    rtdata.get_data("WFI/grism/ramp_grism.asdf")
    rtdata.input = "ramp_grism.asdf"

    # Test CRDS
    step = SaturationStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "saturation")
    ref_file_name = os.path.split(ref_file_path)[-1]

    # Confirm that a WFI saturation reference file name was returned
    assert "roman_wfi_saturation" in ref_file_name

    # Test SaturationStep
    output = "ramp_grism_saturationstep.asdf"
    rtdata.output = output

    args = ["romancal.step.SaturationStep", rtdata.input]
    RomanStep.from_cmdline(args)

    ramp_out = rdm.open(rtdata.output)
    assert ("roman.pixeldq" in ramp_out.to_flat_dict())

    rtdata.get_truth(f"truth/WFI/grism/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)
