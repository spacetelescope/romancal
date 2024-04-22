import asdf
import numpy as np
import pytest
from astropy import units as u
from metrics_logger.decorators import metrics_logger
from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep
from romancal.tweakreg.tweakreg_step import TweakRegStep

from .regtestdata import compare_asdf


@metrics_logger("DMS280", "DMS405", "DMS488")
@pytest.mark.bigdata
def test_tweakreg(rtdata, ignore_asdf_paths, tmp_path):
    # N.B.: uncal file is from simulator
    # ``shifted'' version is created in make_regtestdata.sh; cal file is taken,
    # the wcsinfo is perturbed, and AssignWCS is run to update the WCS with the
    # perturbed information
    orig_uncal = "r0000101001001001001_01101_0001_WFI01_uncal.asdf"
    input_data = "r0000101001001001001_01101_0001_WFI01_shift_cal.asdf"
    output_data = "r0000101001001001001_01101_0001_WFI01_shift_tweakregstep.asdf"
    truth_data = "r0000101001001001001_01101_0001_WFI01_shift_tweakregstep.asdf"

    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.get_data(f"WFI/image/{orig_uncal}")
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = input_data
    rtdata.output = output_data

    # instantiate TweakRegStep (for running and log access)
    step = TweakRegStep()

    args = [
        "romancal.step.TweakRegStep",
        rtdata.input,
        f"--output_file='{rtdata.output}'",
        "--suffix='tweakregstep'",
    ]
    RomanStep.from_cmdline(args)
    tweakreg_out = rdm.open(rtdata.output)

    step.log.info(
        "DMS280 MSG: TweakReg step recorded as complete? :"
        f' {tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"}'
    )
    assert tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"

    step.log.info(
        f"""DMS280 MSG: TweakReg created new attribute with fit results? :\
            {"wcs_fit_results" in tweakreg_out.meta}"""
    )
    assert "wcs_fit_results" in tweakreg_out.meta

    step.log.info(
        f"""DMS280 MSG: TweakReg created new coordinate frame 'v2v3corr'? :\
            {"v2v3corr" in tweakreg_out.meta.wcs.available_frames}"""
    )
    assert "v2v3corr" in tweakreg_out.meta.wcs.available_frames

    diff = compare_asdf(rtdata.output, rtdata.truth, atol=1e-3, **ignore_asdf_paths)
    step.log.info(
        f"DMS280 MSG: Was the proper TweakReg data produced? : {diff.identical}"
    )
    assert diff.identical, diff.report()

    wcstweak = tweakreg_out.meta.wcs
    orig_model_asdf = asdf.open(orig_uncal)
    wcstrue = orig_model_asdf["romanisim"]["wcs"]  # simulated, true WCS
    pts = np.linspace(0, 4000, 30)
    xx, yy = np.meshgrid(pts, pts)
    coordtweak = wcstweak.pixel_to_world(xx, yy)
    coordtrue = wcstrue.pixel_to_world(xx, yy)
    diff = coordtrue.separation(coordtweak).to(u.arcsec).value
    rms = np.sqrt(np.mean(diff**2)) * 1000  # rms difference in mas
    passmsg = "PASS" if rms < 1.3 / np.sqrt(2) else "FAIL"
    step.log.info(
        f"DMS488 MSG: WCS agrees with true WCS to {rms:5.2f} mas, less than "
        f"1.3 / sqrt(2)?  {passmsg}"
    )
    step.log.info(
        f"DMS405 MSG: WCS agrees with true WCS to {rms:5.2f} mas, less than "
        f"5 / sqrt(2)?  {passmsg}"
    )

    assert rms < 1.3 / np.sqrt(2) * 100
