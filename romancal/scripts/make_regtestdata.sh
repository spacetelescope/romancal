#!/bin/bash

# takes one argument: the directory into which to start putting the regtest files.

# input files needed:
# r0000101001001001001_0001_wfi01 - default for most steps
# r0000201001001001001_0001_wfi01 - equivalent for spectroscopic data
# r0000101001001001001_0002_wfi01 - a second resample exposure, only cal step needed
# r0000101001001001001_0003_wfi01 - for ramp fitting; truncated image
# r0000201001001001001_0003_wfi01 - for ramp fitting; truncated spectroscopic
#                                         we need only darkcurrent & ramp fit for these
#
# r00r1601001001001001_0001_wfi01 - special 16 resultant file, imaging, only need cal file
# r10r1601001001001001_0001_wfi01 - special 16 resultant file, spectroscopy, only need cal file


if [ $# -eq 0 ]; then
    echo "Please provide an output directory as a command line argument"
    exit 1
fi

if [ -z "${PATCH_TABLE_PATH}" ]; then
    echo "Please set the PATCH_TABLE_PATH environment variable"
    exit 1
fi

outdir="$1"

logfile="$outdir/make_regtestdata.log"

# stop on an error
# FIXME we can't do this because asn_from_list always returns an error
# see: https://github.com/spacetelescope/romancal/issues/1535
#set -e

# Redirect all output to the logfile
exec > $logfile 2>&1

# set up the directory structure
mkdir -p $outdir/roman-pipeline/dev/WFI/image
mkdir -p $outdir/roman-pipeline/dev/truth/WFI/image
mkdir -p $outdir/roman-pipeline/dev/WFI/grism
mkdir -p $outdir/roman-pipeline/dev/truth/WFI/grism


# most regtests; run the pipeline and save a lot of results.
for fn in r0000101001001001001_0001_wfi01 r0000201001001001001_0001_wfi01
do
    echo "Running pipeline on ${fn}..."
    strun roman_elp ${fn}_uncal.asdf --steps.dq_init.save_results True --steps.saturation.save_results True --steps.linearity.save_results True --steps.dark_current.save_results True --steps.rampfit.save_results True --steps.assign_wcs.save_results True --steps.photom.save_results True --steps.refpix.save_results True --steps.flatfield.save_results True --steps.assign_wcs.save_results True
    [[ ${fn} = r00002* ]] && dirname="grism" || dirname="image"
    cp ${fn}_*.asdf $outdir/roman-pipeline/dev/WFI/$dirname/
    cp ${fn}_*.asdf $outdir/roman-pipeline/dev/truth/WFI/$dirname/

    # uncal file doesn't belong in truth.
    rm $outdir/roman-pipeline/dev/truth/WFI/${dirname}/${fn}_uncal.asdf
done


# truncated files for ramp fit regtests
for fn in r0000101001001001001_0003_wfi01 r0000201001001001001_0003_wfi01
do
    echo "Running pipeline on ${fn}..."
    strun roman_elp ${fn}_uncal.asdf --steps.dark_current.save_results True --steps.rampfit.save_results True
    [[ ${fn} = r00002* ]] && dirname="grism" || dirname="image"
    cp ${fn}_darkcurrent.asdf $outdir/roman-pipeline/dev/WFI/$dirname/
    cp ${fn}_rampfit.asdf $outdir/roman-pipeline/dev/truth/WFI/$dirname/
done
cp r0000101001001001001_0003_wfi01_cal.asdf $outdir/roman-pipeline/dev/WFI/image/


# second imaging exposure
strun roman_elp r0000101001001001001_0002_wfi01_uncal.asdf
cp r0000101001001001001_0002_wfi01_cal.asdf $outdir/roman-pipeline/dev/WFI/image/

# image used in the skycell generation test
strun roman_elp r0000101001001001001_0002_wfi10_uncal.asdf
cp r0000101001001001001_0002_wfi10_cal.asdf $outdir/roman-pipeline/dev/WFI/image/


# CRDS test needs the "usual" r00001..._0001_wfi01 files.
# It also needs a hacked r00001..._0001_wfi01 file, with the time changed.
# this makes the hacked version.
echo "Creating regtest files for CRDS tests..."
basename="r0000101001001001001_0001_wfi01"
python -c "
import asdf
from roman_datamodels import stnode
from astropy.time import Time
basename = '$basename'
f = asdf.open(f'{basename}_uncal.asdf')
f['roman']['meta']['exposure']['start_time'] = Time('2020-01-01T00:00:00', format='isot')
f['roman']['meta']['filename'] = stnode.Filename(f'{basename}_changetime_uncal.asdf')
f.write_to(f'{basename}_changetime_uncal.asdf')"
strun roman_elp ${basename}_changetime_uncal.asdf --steps.assign_wcs.save_results True --steps.flatfield.save_results True
# copy input and truth files into location
cp ${basename}_changetime_assignwcs.asdf $outdir/roman-pipeline/dev/WFI/image
cp ${basename}_changetime_flat.asdf $outdir/roman-pipeline/dev/truth/WFI/image


# need to make a special ALL_SATURATED file for the all saturated test.
echo "Creating regtest files for all saturated tests..."
basename="r0000101001001001001_0001_wfi01"
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
strun romancal.step.DarkCurrentStep r0000101001001001001_0001_wfi01_linearity.asdf --output_file=Test_dark
cp Test_darkcurrent.asdf $outdir/roman-pipeline/dev/truth/WFI/image/


# make a special linearity file with a different suffix
strun romancal.step.LinearityStep r0000101001001001001_0001_wfi01_refpix.asdf --output_file=Test_linearity
cp Test_linearity.asdf $outdir/roman-pipeline/dev/truth/WFI/image/


# we have a test that runs the flat field step directly on an _L1_ spectroscopic
# file and verifies that it gets skipped.
basename="r0000201001001001001_0001_wfi01"
strun romancal.step.FlatFieldStep ${basename}_assignwcs.asdf
cp ${basename}_flat.asdf $outdir/roman-pipeline/dev/truth/WFI/grism/

# make a version of a file with a different pointing
# note this makes the meta['pointing'] data inconsistent with the wcsinfo?
# we haven't updated the filename in these files, but the regtest mechanism
# also doesn't
# update them, and we need to match.
for basename in r0000101001001001001_0001_wfi01 r0000201001001001001_0001_wfi01
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
for basename in r0000101001001001001_0001_wfi01
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


strun roman_elp r0000101001001001001_0004_wfi01_uncal.asdf
strun roman_elp r0000201001001001001_0004_wfi01_uncal.asdf
cp r0000101001001001001_0004_wfi01_uncal.asdf $outdir/roman-pipeline/dev/WFI/image/
cp r0000201001001001001_0004_wfi01_uncal.asdf $outdir/roman-pipeline/dev/WFI/grism/

l3name="r0099101001001001001_F158_visit"
asn_from_list r0000101001001001001_0001_wfi01_cal.asdf r0000101001001001001_0002_wfi01_cal.asdf r0000101001001001001_0003_wfi01_cal.asdf -o L3_regtest_asn.json --product-name $l3name
strun roman_mos L3_regtest_asn.json
cp L3_regtest_asn.json $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# L3 catalog
strun romancal.step.SourceCatalogStep ${l3name}_coadd.asdf
cp ${l3name}_cat.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# L3 on skycell
l3name="r0099101001001001001_r274dp63x31y81_prompt_F158"
asn_from_list r0000101001001001001_0001_wfi01_cal.asdf r0000101001001001001_0002_wfi01_cal.asdf r0000101001001001001_0003_wfi01_cal.asdf -o L3_mosaic_asn.json --product-name $l3name --target r274dp63x31y81
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
strun romancal.step.ResampleStep L3_mosaic_asn.json --rotation=0 --output_file=mosaic.asdf
cp mosaic_resamplestep.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# L3 catalog
strun romancal.step.SourceCatalogStep ${l3name}_coadd.asdf
cp ${l3name}_cat.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# L2 catalog
strun romancal.step.SourceCatalogStep r0000101001001001001_0001_wfi01_cal.asdf
cp r0000101001001001001_0001_wfi01_cat.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp r0000101001001001001_0001_wfi01_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# L3 skycell catalog
strun romancal.step.SourceCatalogStep ${l3name}_coadd.asdf
cp ${l3name}_cat.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# multiband catalog
asn_from_list --product-name=${l3name}_mbcat ${l3name}_coadd.asdf -o L3_skycell_mbcat_asn.json
strun romancal.step.MultibandCatalogStep L3_skycell_mbcat_asn.json --deblend True
cp L3_skycell_mbcat_asn.json $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_mbcat_cat.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
cp ${l3name}_mbcat_segm.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# tests passing suffix to the pipeline
strun roman_elp r0000101001001001001_0001_wfi01_uncal.asdf --steps.tweakreg.skip=True --suffix=star
cp r0000101001001001001_0001_wfi01_star.asdf $outdir/roman-pipeline/dev/truth/WFI/image/

# 2nd L3 on skycell
l3name="r0099101001001001001_0001_r274dp63x31y81_prompt_F158"
asn_from_list r0000101001001001001_0001_wfi01_cal.asdf -o L3_mosaic_0001_asn.json --product-name $l3name --target r274dp63x31y81
# The pipeline will silently do nothing and not return an error exit code if the output
# file already exists.
# see: https://github.com/spacetelescope/romancal/issues/1544
# To work around this remove the expected output file
if [ -f "${l3name}_coadd.asdf" ]; then
    rm "${l3name}_coadd.asdf"
fi
strun roman_mos L3_mosaic_0001_asn.json
cp L3_mosaic_0001_asn.json $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/WFI/image/
cp ${l3name}_coadd.asdf $outdir/roman-pipeline/dev/truth/WFI/image/
