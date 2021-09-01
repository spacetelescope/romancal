import os
import pytest

from romancal.stpipe import RomanStep
from romancal.step import FlatFieldStep
import roman_datamodels as rdm
from .regtestdata import compare_asdf

@pytest.mark.bigdata
def test_flat_field_image_step(rtdata, ignore_asdf_paths):

    rtdata.get_data("WFI/image/l2_0004_rate.asdf")
    rtdata.input = "l2_0004_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "flat")
    ref_file_name = os.path.split(ref_file_path)[-1]
    assert "roman_wfi_flat" in ref_file_name


    # Test FlatFieldStep
    output = "l2_0004_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)

@pytest.mark.skip(reason="There are no grism flats.")
@pytest.mark.bigdata
def test_flat_field_grism_step(rtdata, ignore_asdf_paths):

    rtdata.get_data("WFI/grism/l2_0004_grism_rate.asdf")
    rtdata.input = "l2_0004_grism_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = rdm.open(rtdata.input)
    ref_file_path = step.get_reference_file(model, "flat")
    ref_file_name = os.path.split(ref_file_path)[-1]
    assert "roman_wfi_flat" in ref_file_name

    # Test FlatFieldStep
    output = "l2_0004_grism_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/grism/{output}")
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)

@pytest.mark.bigdata
def test_flat_field_crds_match_image_step(rtdata, ignore_asdf_paths):
    # DMS79 Test: Testing that different datetimes pull different
    # flat files and successfully make level 2 output

    # First file
    rtdata.get_data("WFI/image/l2_0004_rate.asdf")
    rtdata.input = "l2_0004_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = rdm.open(rtdata.input)
    step.log.info('DMS79 MSG: Testing retrieval of best ref file, '
                  'Success is flat file with correct use after date')


    step.log.info(f'DMS79 MSG: First data file: {rtdata.input.rsplit("/", 1)[1]}')
    step.log.info(f'DMS79 MSG: Observation date: {model.meta.observation.start_time}')

    ref_file_path = step.get_reference_file(model, "flat")
    step.log.info(f'DMS79 MSG: CRDS matched flat file: {ref_file_path.rsplit("/", 1)[1]}')
    flat = rdm.open(ref_file_path)
    step.log.info(f'DMS79 MSG: flat file UseAfter date: {flat.meta.useafter}')
    step.log.info(f'DMS79 MSG: UseAfter date before observation date? : '
                  f'{(flat.meta.useafter < model.meta.observation.start_time)}')

    # Test FlatFieldStep
    output = "l2_0004_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    step.log.info('DMS79 MSG: Running flat fielding step. The first ERROR is expected, '
                  'due to extra CRDS parameters not having been implemented yet.')
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")

    step.log.info(f'DMS79 MSG: Was proper flat fielded Level 2 data produced? : '
                  f'{(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}')
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)

    # This test requires a second file, in order to meet the DMS79 requirement.
    # The test will show that two files with different observation dates match
    #  to separate flat files in CRDS.

    # Second file
    rtdata.get_data("WFI/image/l2_0004b_rate.asdf")
    rtdata.input = "l2_0004b_rate.asdf"

    # Test CRDS
    step = FlatFieldStep()
    model = rdm.open(rtdata.input)

    step.log.info(f'DMS79 MSG: Second data file: {rtdata.input.rsplit("/", 1)[1]}')
    step.log.info(f'DMS79 MSG: Observation date: {model.meta.observation.start_time}')

    ref_file_path_b = step.get_reference_file(model, "flat")
    step.log.info(f'DMS79 MSG: CRDS matched flat file: {ref_file_path_b.rsplit("/", 1)[1]}')
    flat = rdm.open(ref_file_path_b)
    step.log.info(f'DMS79 MSG: flat file UseAfter date: {flat.meta.useafter}')
    step.log.info(f'DMS79 MSG: UseAfter date before observation date? : '
                  f'{(flat.meta.useafter < model.meta.observation.start_time)}')

    # Test FlatFieldStep
    output = "l2_0004b_rate_flatfieldstep.asdf"
    rtdata.output = output
    args = ["romancal.step.FlatFieldStep", rtdata.input]
    step.log.info('DMS79 MSG: Running flat fielding step. The first ERROR is expected, '
                  'due to extra CRDS parameters not having been implemented yet.')
    RomanStep.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    step.log.info(f'DMS79 MSG: Was proper flat fielded Level 2 data produced? : '
                  f'{(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}')
    assert (compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)

    # Test differing flat matches
    step.log.info(f'DMS79 MSG REQUIRED TEST: Are the two data files matched to different '
                  f'flat files? : '
                  f'{("/".join(ref_file_path.rsplit("/", 3)[1:])) != ("/".join(ref_file_path_b.rsplit("/", 3)[1:]))}')
    assert ("/".join(ref_file_path.rsplit("/", 1)[1:])) != \
           ("/".join(ref_file_path_b.rsplit("/", 1)[1:]))
