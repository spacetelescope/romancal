#!/bin/bash

# takes one argument: the directory into which to start putting the regtest files.

# input files needed:
# r0000101001001001001_0001_wfi*_f158_uncal.asdf - L1 files for each detector


if [ $# -eq 0 ]; then
    echo "Please provide an output directory as a command line argument"
    exit 1
fi

outdir="$1"

logfile="$outdir/make_regtestdata.log"

# stop on an error
set -e
mkdir -p $outdir/L1
mkdir -p $outdir/MOC

# set up the directory structure

# Redirect all output to the logfile and the terminal
exec > >(tee $logfile) 2>&1

# run elp for 1 uncal for each wfi
L1_FNS=`ls r0000101001001001001_0001_wfi*_f158_uncal.asdf`
for L1_FN in $L1_FNS; do
    cp $L1_FN $outdir/L1/
    strun roman_elp $L1_FN
done

# copy cal files to output directory
L2_FNS=`ls r0000101001001001001_0001_wfi*_f158_cal.asdf`
for L2_FN in $L2_FNS; do
    cp $L2_FN $outdir/MOC/
done
