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

# Redirect all output to the logfile and the terminal
exec > >(tee $logfile) 2>&1

# set up the directory structure
#mkdir -p $outdir/MOC

# run elp for 1 uncal for each wfi
L1_FNS=`ls r0000101001001001001_0001_wfi*_f158_uncal.asdf`
for L1_FN in $L1_FNS; do
    strun roman_elp $L1_FN
done

# copy cal files to MOC directory
#L2_FNS=`ls r0000101001001001001_0001_wfi*_f158_cal.asdf`
