import os
import asdf
import numpy as np
import pytest
from astropy import units as u
from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf


@pytest.mark.bigdata
def test_tweakreg(
    rtdata, ignore_asdf_paths, tmp_path, resource_tracker, request, dms_logger
):
    # N.B.: uncal file is from simulator
    # ``shifted'' version is created in make_regtestdata.sh; cal file is taken,
    # the wcsinfo is perturbed, and AssignWCS is run to update the WCS with the
    # perturbed information
    orig_uncal = "r0000101001001001001_0001_wfi01_f158_uncal.asdf"
    orig_catfile = "r0000101001001001001_0001_wfi01_f158_cat.parquet"
    input_data = "r0000101001001001001_0001_wfi01_f158_shift_cal.asdf"
    output_data = "r0000101001001001001_0001_wfi01_f158_shift_tweakregstep.asdf"
    truth_data = "r0000101001001001001_0001_wfi01_f158_shift_tweakregstep.asdf"

    rtdata.get_data(f"WFI/image/{input_data}")
    rtdata.get_data(f"WFI/image/{orig_catfile}")
    rtdata.get_data(f"WFI/image/{orig_uncal}")
    rtdata.get_truth(f"truth/WFI/image/{truth_data}")

    rtdata.input = input_data
    rtdata.output = output_data

    args = [
        "romancal.step.TweakRegStep",
        rtdata.input,
        f"--output_file='{rtdata.output}'",
        "--suffix='tweakregstep'",
    ]
    with resource_tracker.track(log=request):
        RomanStep.from_cmdline(args)

    tweakreg_out = rdm.open(rtdata.output)

    dms_logger.info(
        "DMS280 MSG: TweakReg step recorded as complete? :"
        f" {tweakreg_out.meta.cal_step.tweakreg == 'COMPLETE'}"
    )
    assert tweakreg_out.meta.cal_step.tweakreg == "COMPLETE"

    dms_logger.info(
        f"""DMS280 MSG: TweakReg created new attribute with fit results? :\
            {"wcs_fit_results" in tweakreg_out.meta}"""
    )
    assert "wcs_fit_results" in tweakreg_out.meta

    dms_logger.info(
        f"""DMS280 MSG: TweakReg created new coordinate frame 'v2v3corr'? :\
            {"v2v3corr" in tweakreg_out.meta.wcs.available_frames}"""
    )
    assert "v2v3corr" in tweakreg_out.meta.wcs.available_frames

    wcstweak = tweakreg_out.meta.wcs
    with asdf.open(orig_uncal) as orig_model_asdf:
        wcstrue = orig_model_asdf["romanisim"]["wcs"]  # simulated, true WCS
        pts = np.linspace(0, 4000, 30)
        xx, yy = np.meshgrid(pts, pts)
        coordtweak = wcstweak.pixel_to_world(xx, yy)
        coordtrue = wcstrue.pixel_to_world(xx, yy)
    diff = coordtrue.separation(coordtweak).to(u.arcsec).value
    rms = np.sqrt(np.mean(diff**2)) * 1000  # rms difference in mas
    passmsg = "PASS" if rms < 1.3 / np.sqrt(2) else "FAIL"
    dms_logger.info(
        f"DMS488 MSG: WCS agrees with true WCS to {rms:5.2f} mas, less than "
        f"1.3 / sqrt(2)?  {passmsg}"
    )
    dms_logger.info(
        f"DMS405 MSG: WCS agrees with true WCS to {rms:5.2f} mas, less than "
        f"5 / sqrt(2)?  {passmsg}"
    )

    assert rms < 1.3 / np.sqrt(2)

    #check that tweakreg has been run
    if tweakreg_out.meta.cal_step.tweakreg == "COMPLETE":
        # Find the reference catalog used by tweakreg
        for s in tweakreg_out.meta.cal_logs:
            if 'abs_refcat' in s:
                str = s

            str1 = str[str.rfind('abs_refcat'):]
            refcat_name = str1[:str1.find('\n')]

        passmsg = "PASS" if "GAIA" in refcat_name else "FAIL"
        dms_logger.info(f"DMS549 MSG: {passmsg}, {refcat_name} used to align data to"
                        f"the Gaia astrometric reference frame.")
    
    # check if the Mean Absolute Error is less that 5 milliarcsec (10.0E-3 * pixelscale)
    assert  tweakreg_out.meta.wcs_fit_results.mae * 0.1 < 5.0e-3 
    passmsg = "PASS" if  tweakreg_out.meta.wcs_fit_results.mae * 0.1 < 5.0e-3 else "FAIL"
    dms_logger.info(f"DMS406 {passmsg}, the Mean Absolute Error is less that 5 milliarcsec.")
    # check that the tweakreg step is marked complete
    passmsg = "PASS" if  tweakreg_out.meta.cal_step.tweakreg == "COMPLETE" else "FAIL"
    dms_logger.info(f"DMS406 {passmsg}, the Tweakreg step is compete.")
    wcs_filename = output_data.rsplit("_", 1)[0] + "_wcs.asdf"
    passmsg = "PASS" if os.path.isfile(wcs_filename) else "FAIL"
    dms_logger.info(f"DMS406 {passmsg}, the output wcs file exists.")
    diff = compare_asdf(rtdata.output, rtdata.truth, atol=1e-3, **ignore_asdf_paths)
    dms_logger.info(
        f"DMS280 MSG: Was the proper TweakReg data produced? : {diff.identical}"
    )
    assert diff.identical, diff.report()
