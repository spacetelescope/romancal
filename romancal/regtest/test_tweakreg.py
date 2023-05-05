from roman_datamodels import datamodels as rdm
from romancal.tweakreg.tweakreg_step import TweakRegStep
from romancal.stpipe import RomanStep
from roman_datamodels import maker_utils
import os
import csv
import asdf
import pytest
from .regtestdata import compare_asdf


def load_base_image_wcs(input_dm):
    # data from romanisim simulated image 'l1-270-66-gaia-2016-sca7_cal.asdf'
    # command used to create simulated image:
    # romanisim-make-image --catalog ./rsim_cat_F158.ecsv --radec 270.0 66.0 --bandpass F158 --sca 7
    # --usecrds --webbpsf --date 2026 1 1 --level 1 l1-270-66-gaia-2016-sca7.asdf

    # update meta.wcsinfo
    input_dm.meta.wcsinfo.v2_ref = 0.42955128282521254
    input_dm.meta.wcsinfo.v3_ref = -0.2479976768255853
    input_dm.meta.wcsinfo.vparity = -1
    input_dm.meta.wcsinfo.ra_ref = 270.0
    input_dm.meta.wcsinfo.dec_ref = 66.0
    input_dm.meta.wcsinfo.roll_ref = 60

    # read saved WCS object
    full_path_to_wcs_file = os.path.join(
        os.path.dirname(__file__), "data/base_image_wcs.asdf"
    )
    asdf_file = asdf.open(full_path_to_wcs_file)
    wcs = asdf_file.tree["wcs"]

    return wcs


def create_base_image_source_catalog(tmpdir, output_filename):
    # create a temp CSV file to be used as source catalog
    header = ["x", "y"]
    data = [
        [100.00, 50.00],
        [200.00, 100.00],
        [400.00, 200.00],
    ]
    output = os.path.join(tmpdir, output_filename)
    with open(output, "w", encoding="UTF8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)


def create_base_image_and_catalog(tmpdir):
    # data from romanisim simulated image 'l1-270-66-gaia-2016-sca7_cal.asdf'
    l2 = maker_utils.mk_level2_image(shape=(4088, 4088))

    # update wcsinfo and add WCS from a simulated image
    l2.meta["wcs"] = load_base_image_wcs(l2)

    # add source detection catalog name
    tweakreg_catalog_filename = "base_image_sources.csv"
    create_base_image_source_catalog(tmpdir, tweakreg_catalog_filename)
    l2.meta["tweakreg_catalog"] = os.path.join(
        tmpdir, tweakreg_catalog_filename
    )

    l2im = rdm.ImageModel(l2)
    return l2im


@pytest.mark.bigdata
def test_is_wcs_correction_small(rtdata, ignore_asdf_paths):
    input_data = "r0000401001001001001_01101_0001_WFI01_cal_tweakreg_input.asdf"
    output_data = (
        "r0000401001001001001_01101_0001_WFI01_cal_tweakreg_output.asdf"
    )
    truth_data = "r0000401001001001001_01101_0001_WFI01_cal_tweakreg_truth.asdf"

    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = input_data
    rtdata.output = output_data

    step = TweakRegStep()

    args = ["romancal.step.TweakRegStep", rtdata.input]
    RomanStep.from_cmdline(args)
    tweakreg_out = rdm.open(rtdata.output)

    step.log.info(
        "DMSXXX MSG: TweakReg step recorded as complete? :"
        f' {tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"}'
    )
    assert tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"

    step.log.info(
        "DMSXXX MSG: Was the proper TweakReg data produced?"
        f" : {(compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None)}"
    )
    assert (
        compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths) is None
    )
