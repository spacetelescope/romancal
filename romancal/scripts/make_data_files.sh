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

# set up the directory structure
mkdir -p $outdir/L1
mkdir -p $outdir/MOC
mkdir -p $outdir/data-workshop

# Redirect all output to the logfile and the terminal
exec > >(tee $logfile) 2>&1

# run elp for 1 uncal image for each wfi
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

# data-workshop files
# L1 -> L2
WORKSHOP_UNCALS="
    r0000101001001001001_0001_wfi01_f158_uncal.asdf
    r0000201001001001001_0001_wfi01_grism_uncal.asdf
    r0000301001001001001_0001_wfi01_prism_uncal.asdf
"
for WORKSHOP_UNCAL in $WORKSHOP_UNCALS; do
    if [ $WORKSHOP_UNCAL =~ "_f158_"]; then
        # image file, already processed through elp to a cal file above
        cp "${WORKSHOP_UNCAL//_uncal/_cal}" $outdir/data-workshop/
        cp "${WORKSHOP_UNCAL//_uncal/_wcs}" $outdir/data-workshop/
        cp "${WORKSHOP_UNCAL//_uncal/_segm}" $outdir/data-workshop/
        cp "${WORKSHOP_UNCAL//_uncal.asdf/_cat.parquet}" $outdir/data-workshop/

        # prepare associations
        skycell_asn "${WORKSHOP_UNCAL//_uncal/_cal}" -o r00001 --product-type visit
    else
        # grism/prism file, only run through elp
        strun roman_elp $WORKSHOP_UNCAL
        cp "${WORKSHOP_UNCAL//_uncal/_cal}" $outdir/data-workshop/
        cp "${WORKSHOP_UNCAL//_uncal/_wcs}" $outdir/data-workshop/
    fi
done

# L2 -> L3
ASN_FNS=`ls *_asn.json`
for ASN_FN in $ASN_FNS; do
    strun roman_mos $ASN_FN
    cp $ASN_FN $outdir/data-workshop/
    cp "${WORKSHOP_UNCAL//_asn.json/_coadd.asdf}" $outdir/data-workshop/
    cp "${WORKSHOP_UNCAL//_asn.json/_cat.parquet}" $outdir/data-workshop/
    cp "${WORKSHOP_UNCAL//_asn.json/_segm.asdf}" $outdir/data-workshop/
done
