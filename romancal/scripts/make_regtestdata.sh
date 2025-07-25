#!/bin/bash

# takes one argument: the directory into which to start putting the regtest files.

# input files needed:
# r0000101001001001001_0001_wfi01_f158  - default for most steps
# r0000201001001001001_0001_wfi01_grism - equivalent for spectroscopic data
# r0000101001001001001_0002_wfi01_f158  - a second resample exposure, only cal step needed
# r0000101001001001001_0003_wfi01_f158  - for ramp fitting; truncated image
# r0000201001001001001_0003_wfi01_grism - for ramp fitting; truncated spectroscopic
#                                         we need only darkcurrent & ramp fit for these
# r0000101001001001001_0004_wfi01_f158  - 16 resultant imaging file
# r0000201001001001001_0004_wfi01_grism - 16 resultant spectroscopy file

# this script also downloads two files from artifactory for TVAC tests; these
# do not come from romanisim.

if [ $# -eq 0 ]; then
    echo "Please provide an output directory as a command line argument"
    exit 1
fi

outdir="$1"

logfile="$outdir/make_regtestdata.log"

# stop on an error
set -e

# Redirect all output to the logfile and the terminal
exec > >(tee $logfile) 2>&1

# set up the directory structure
mkdir -p $outdir/roman-pipeline/dev/WFI/image
mkdir -p $outdir/roman-pipeline/dev/truth/WFI/image
mkdir -p $outdir/roman-pipeline/dev/WFI/grism
mkdir -p $outdir/roman-pipeline/dev/truth/WFI/grism
mkdir -p $outdir/roman-pipeline/dev/references


# most regtests; run the pipeline and save a lot of results.
for fn in r0000101001001001001_0001_wfi01_f158 r0000201001001001001_0001_wfi01_grism
do
    echo "Running pipeline on ${fn}..."
    strun roman_elp ${fn}_uncal.asdf --steps.dq_init.save_results True --steps.saturation.save_results True --steps.linearity.save_results True --steps.dark_current.save_results True --steps.rampfit.save_results True --steps.assign_wcs.save_results True --steps.photom.save_results True --steps.refpix.save_results True --steps.flatfield.save_results True --steps.assign_wcs.save_results True
    [[ ${fn} = r00002* ]] && dirname="grism" || dirname="image"
    cp ${fn}_*.asdf $outdir/roman-pipeline/dev/WFI/$dirname/
    cp ${fn}_*.asdf $outdir/roman-pipeline/dev/truth/WFI/$dirname/

    # uncal file doesn't belong in truth.
    rm $outdir/roman-pipeline/dev/truth/WFI/${dirname}/${fn}_uncal.asdf
done

# L2 catalog
cp r0000101001001001001_0001_wfi01_f158_cat.parquet $outdir/roman-pipeline/dev/truth/WFI/image/
cp r0000101001001001001_0001_wfi01_f158_cat.parquet $outdir/roman-pipeline/dev/WFI/image/

# truncated files for ramp fit regtests
for fn in r0000101001001001001_0003_wfi01_f158 r0000201001001001001_0003_wfi01_grism
do
    echo "Running pipeline on ${fn}..."
    strun roman_elp ${fn}_uncal.asdf --steps.dark_current.save_results True --steps.rampfit.save_results True
    [[ ${fn} = r00002* ]] && dirname="grism" || dirname="image"
    cp ${fn}_darkcurrent.asdf $outdir/roman-pipeline/dev/WFI/$dirname/
    cp ${fn}_rampfit.asdf $outdir/roman-pipeline/dev/truth/WFI/$dirname/
done
cp r0000101001001001001_0003_wfi01_f158_cal.asdf $outdir/roman-pipeline/dev/WFI/image/


# second imaging exposure
strun roman_elp r0000101001001001001_0002_wfi01_f158_uncal.asdf
cp r0000101001001001001_0002_wfi01_f158_cal.asdf $outdir/roman-pipeline/dev/WFI/image/

# image used in the skycell generation test
strun roman_elp r0000101001001001001_0002_wfi10_f158_uncal.asdf
cp r0000101001001001001_0002_wfi10_f158_cal.asdf $outdir/roman-pipeline/dev/WFI/image/

# need to make a special ALL_SATURATED file for the all saturated test.
echo "Creating regtest files for all saturated tests..."
basename="r0000101001001001001_0001_wfi01_f158"
python -c "
import asdf
from roman_datamodels import stnode
basename = '$basename'
f = asdf.open(f'{basename}_uncal.asdf')
data = f['roman']['data'].copy()
data[...] = 65535
f['roman']['data'] = data
f['roman']['meta']['filename'] = stnode.Filename(f'{basename}_ALL_SATURATED_uncal.asdf')
f.write_to(f'{basename}_ALL_SATURATED_uncal.asdf')"
strun roman_elp ${basename}_ALL_SATURATED_uncal.asdf
cp ${basename}_ALL_SATURATED_uncal.asdf $outdir/roman-pipeline/dev/WFI/image/
cp ${basename}_ALL_SATURATED_cal.asdf $outdir/roman-pipeline/dev/truth/WFI/image/


# make a special file dark file with a different name
strun romancal.step.DarkCurrentStep r0000101001001001001_0001_wfi01_f158_rampfit.asdf --output_file=Test_dark
cp Test_darkcurrent.asdf $outdir/roman-pipeline/dev/truth/WFI/image/


# make a special linearity file with a different suffix
strun romancal.step.LinearityStep r0000101001001001001_0001_wfi01_f158_refpix.asdf --output_file=Test_linearity
cp Test_linearity.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# test likelihood-based ramp fitting step
strun romancal.step.RampFitStep r0000101001001001001_0001_wfi01_f158_linearity.asdf --algorithm=likely --output_file=r0000101001001001001_0001_wfi01_f158_like_rampfit.asdf
cp r0000101001001001001_0001_wfi01_f158_like_rampfit.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# We have a test that runs the flat field step directly on an _L1_ spectroscopic
# file and verifies that it gets skipped.
basename="r0000201001001001001_0001_wfi01_grism"
strun romancal.step.FlatFieldStep ${basename}_assignwcs.asdf
cp ${basename}_flat.asdf $outdir/roman-pipeline/dev/truth/WFI/grism/

# make a version of a file with a different pointing
# note this makes the meta['pointing'] data inconsistent with the wcsinfo?
# we haven't updated the filename in these files, but the regtest mechanism
# also doesn't
# update them, and we need to match.
for basename in r0000101001001001001_0001_wfi01_f158 r0000201001001001001_0001_wfi01_grism
do
    python -c "
import asdf
import roman_datamodels as rdm
from romancal.assign_wcs.assign_wcs_step import AssignWcsStep
model = rdm.open('${basename}_cal.asdf', lazy_load=False)
delta = [10.0, 10.0]
wcsinfo = model.meta.wcsinfo
if wcsinfo.ra_ref >= 350.0:
    delta[0] *= -1.0
if wcsinfo.dec_ref >= 80.0:
    delta[1] *= -1.0
wcsinfo.ra_ref += delta[0]
wcsinfo.dec_ref += delta[1]
model = AssignWcsStep.call(model)
model.to_asdf(f'${basename}_cal_repoint.asdf')"
    [[ ${basename} = r00002* ]] && dirname="grism" || dirname="image"
    cp ${basename}_cal_repoint.asdf $outdir/roman-pipeline/dev/truth/WFI/$dirname/
done

# Test tweakreg with repointed file, only shifted by 1"
for basename in r0000101001001001001_0001_wfi01_f158
do
    python -c "
import asdf
import roman_datamodels as rdm
from roman_datamodels import stnode
from romancal.assign_wcs.assign_wcs_step import AssignWcsStep
model = rdm.open('${basename}_cal.asdf', lazy_load=False)
model.meta.filename = stnode.Filename(f'${basename}_shift_cal.asdf')
delta = [1 / 3600., 1 / 3600.]
wcsinfo = model.meta.wcsinfo
wcsinfo.ra_ref += delta[0]
wcsinfo.dec_ref += delta[1]
model = AssignWcsStep.call(model)
model.to_asdf(f'${basename}_shift_cal.asdf')"
    strun romancal.step.TweakRegStep ${basename}_shift_cal.asdf
    cp ${basename}_shift_cal.asdf $outdir/roman-pipeline/dev/WFI/image/
    cp ${basename}_shift_tweakregstep.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
done


strun roman_elp r0000101001001001001_0004_wfi01_f158_uncal.asdf
strun roman_elp r0000201001001001001_0004_wfi01_grism_uncal.asdf
cp r0000101001001001001_0004_wfi01_f158_uncal.asdf $outdir/roman-pipeline/dev/WFI/image/
cp r0000201001001001001_0004_wfi01_grism_uncal.asdf $outdir/roman-pipeline/dev/WFI/grism/

# tests passing suffix to the pipeline
strun roman_elp r0000101001001001001_0001_wfi01_f158_uncal.asdf --steps.tweakreg.skip=True --suffix=star
cp r0000101001001001001_0001_wfi01_f158_star.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

l3name="r0000101001001001001_f158"
asn_from_list r0000101001001001001_0001_wfi01_f158_cal.asdf r0000101001001001001_0002_wfi01_f158_cal.asdf r0000101001001001001_0003_wfi01_f158_cal.asdf -o L3_regtest_asn.json --product-name $l3name
strun roman_mos L3_regtest_asn.json
cp L3_regtest_asn.json $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_cat.parquet $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# L3 on skycell
l3name="r00001_p_v01001001001001_270p65x49y70_f158"
asn_from_list r0000101001001001001_0001_wfi01_f158_cal.asdf r0000101001001001001_0002_wfi01_f158_cal.asdf r0000101001001001001_0003_wfi01_f158_cal.asdf -o L3_mosaic_asn.json --product-name $l3name --target 270p65x49y70
# The pipeline will silently do nothing and not return an error exit code if the output
# file already exists.
# see: https://github.com/spacetelescope/romancal/issues/1544
# To work around this remove the expected output file
if [ -f "${l3name}_coadd.asdf" ]; then
    rm "${l3name}_coadd.asdf"
fi
strun roman_mos L3_mosaic_asn.json
cp L3_mosaic_asn.json $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_cat.parquet $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
# also copy outside of truth for input to forced photometry tests
cp ${l3name}_cat.parquet $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_segm.asdf $outdir/roman-pipeline/dev/WFI/image/

strun romancal.step.ResampleStep L3_mosaic_asn.json --resample_on_skycell=False --rotation=0 --output_file=mosaic.asdf
cp mosaic_resamplestep.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# multiband catalog
asn_from_list --product-name=${l3name}_mbcat ${l3name}_coadd.asdf -o L3_skycell_mbcat_asn.json
strun romancal.step.MultibandCatalogStep L3_skycell_mbcat_asn.json --deblend True
cp L3_skycell_mbcat_asn.json $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_mbcat_cat.parquet $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_mbcat_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# 2nd L3 on skycell
l3name="r00001_p_e01001001001001_0001_270p65x49y70_f158"
asn_from_list r0000101001001001001_0001_wfi01_f158_cal.asdf -o L3_mosaic_0001_asn.json --product-name $l3name --target 270p65x49y70
# The pipeline will silently do nothing and not return an error exit code if the output
# file already exists.
# see: https://github.com/spacetelescope/romancal/issues/1544
# To work around this remove the expected output file
if [ -f "${l3name}_coadd.asdf" ]; then
    rm "${l3name}_coadd.asdf"
fi
strun roman_mos L3_mosaic_0001_asn.json
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/WFI/image/

# forced photometry on shallow skycell from deep skycell
strun romancal.step.SourceCatalogStep ${l3name}_coadd.asdf --forced_segmentation r00001_p_v01001001001001_270p65x49y70_f158_segm.asdf --output_file ${l3name}_force_cat.parquet
cp ${l3name}_force_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_force_cat.parquet $outdir/roman-pipeline/dev/truth/WFI/image/

jf rt dl roman-pipeline/dev/WFI/image/TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_uncal.asdf --flat
jf rt dl roman-pipeline/dev/references/dark_ma510.asdf --flat
strun roman_elp TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_uncal.asdf --steps.tweakreg.skip=true --steps.source_catalog.skip=true --steps.dq_init.save=true --steps.dark_current.override_dark=dark_ma510.asdf
cp TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_uncal.asdf $outdir/roman-pipeline/dev/WFI/image/
cp TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_cal.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp TVAC2_NOMOPS_WFIFLA_20240419194120_WFI01_dqinit.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp dark_ma510.asdf $outdir/roman-pipeline/dev/references/
