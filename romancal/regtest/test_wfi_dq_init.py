import os
import pytest

from romancal.stpipe import RomanStep
from romancal.step import DQInitStep
import roman_datamodels as rdm
from .regtestdata import compare_asdf

@pytest.mark.bigdata
def test_dq_init_image_step(rtdata, ignore_asdf_paths):
    # DMS25 Test: Testing retrieval of best ref file for image data,
    # and creation of a ramp file with CRDS selected mask file applied.

    rtdata.get_data("WFI/image/l1_0005_science_raw.asdf")
    rtdata.input = "l1_0005_science_raw.asdf"

    # Test CRDS
    step = DQInitStep()
    model = rdm.open(rtdata.input)
    step.log.info('DMS25 MSG: Testing retrieval of best ref file for image data, '
                  'Success is creation of a ramp file with CRDS selected mask file applied.')

    step.log.info(f'DMS25 MSG: First data file: {rtdata.input.rsplit("/", 1)[1]}')
    ref_file_path = step.get_reference_file(model, "mask")
    step.log.info(f'DMS25 MSG: CRDS matched mask file: {ref_file_path.rsplit("/", 1)[1]}')
    ref_file_name = os.path.split(ref_file_path)[-1]

    assert "roman_wfi_mask" in ref_file_name

    # Test DQInitStep
    output = "l1_0005_science_raw_dqinitstep.asdf"
    rtdata.output = output
    args = ["romancal.step.DQInitStep", rtdata.input]
    step.log.info('DMS25 MSG: Running data quality initialization step. The first ERROR is '
                  'expected, due to extra CRDS parameters not having been implemented yet.')
    RomanStep.from_cmdline(args)
    ramp_out = rdm.open(rtdata.output)
    step.log.info(f'DMS25 MSG: Does ramp data contain pixeldq from mask file? : '
                  f'{("roman.pixeldq" in ramp_out.to_flat_dict())}')
    assert ("roman.pixeldq" in ramp_out.to_flat_dict())

    rtdata.get_truth(f"truth/WFI/image/{output}")
    step.log.info(f'DMS25 MSG: Was proper data quality initialized ramp data produced? : '
                  f'{(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}')
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)




@pytest.mark.bigdata
def test_dq_init_grism_step(rtdata, ignore_asdf_paths):
    # DMS26 Test: Testing retrieval of best ref file for grism data,
    # and creation of a ramp file with CRDS selected mask file applied.

    rtdata.get_data("WFI/grism/l1_0005_grism_science_raw.asdf")
    rtdata.input = "l1_0005_grism_science_raw.asdf"

    # Test CRDS
    step = DQInitStep()
    model = rdm.open(rtdata.input)
    step.log.info('DMS26 MSG: Testing retrieval of best ref file for grism data, '
                  'Success is creation of a ramp file with CRDS selected mask file applied.')

    step.log.info(f'DMS26 MSG: First data file: {rtdata.input.rsplit("/", 1)[1]}')
    ref_file_path = step.get_reference_file(model, "mask")
    step.log.info(f'DMS26 MSG: CRDS matched mask file: {ref_file_path.rsplit("/", 1)[1]}')
    ref_file_name = os.path.split(ref_file_path)[-1]

    assert "roman_wfi_mask" in ref_file_name

    # Test DQInitStep
    output = "l1_0005_grism_science_raw_dqinitstep.asdf"
    rtdata.output = output
    args = ["romancal.step.DQInitStep", rtdata.input]
    step.log.info('DMS26 MSG: Running data quality initialization step. The first ERROR is '
                  'expected, due to extra CRDS parameters not having been implemented yet.')
    RomanStep.from_cmdline(args)
    ramp_out = rdm.open(rtdata.output)
    step.log.info(f'DMS26 MSG: Does ramp data contain pixeldq from mask file? : '
                  f'{("roman.pixeldq" in ramp_out.to_flat_dict())}')
    assert ("roman.pixeldq" in ramp_out.to_flat_dict())

    rtdata.get_truth(f"truth/WFI/grism/{output}")
    step.log.info(f'DMS26 MSG: Was proper data quality initialized ramp data produced? : '
                  f'{(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}')
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)
