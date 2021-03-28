import os
import pytest
from asdf.commands import diff as asdf_diff

from romancal.stpipe import RomanStep
from romancal.step import FlatFieldStep
from romancal import datamodels


@pytest.mark.bigdata
def test_flat_field_step(rtdata):

    rtdata.get_data("WFI/image/l2_0001_rate.asdf")
    rtdata.input = "l2_0001_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = datamodels.ImageModel(rtdata.input)
    ref_file_path = step.get_reference_file(model, "flat")
    ref_file_name = os.path.split(ref_file_path)[-1]
    assert "roman_wfi_flat" in ref_file_name

    # Test FlatFieldStep
    output = "l2_0001_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    compare(rtdata.output, rtdata.truth, **kwargs)


def compare(result, truth, **kwargs):
    f = StringIO()
    asdf_diff([rtdata.output, rtdata.truth], minimal=False,
               iostream=StringIO(), **wkargs)
    if f.getavlue():
        f.get_value()
