Description
============
The Data Quality (DQ) initialization step in the calibration pipeline performs
two functions: Population of the DQ mask and the conversion of the science raw
data, the "_uncal" file, to the ramp model by the majority of steps in the
exposure level pipeline.

Data Quality Initialization
---------------------------

The Data Quality (DQ) initialization step in the calibration pipeline
populates the DQ mask for the input dataset. Flag values from the
appropriate static mask ("MASK") reference file in CRDS are copied into the
"PIXELDQ" array of the input dataset, because it is assumed that flags in the
mask reference file pertain to problem conditions that affect all groups for
a given pixel.

The actual process consists of the following steps:

 - Determine what MASK reference file to use via the interface to the bestref
   utility in CRDS.

 - Copy the input product into a RampModel (if it isn't already) for processing
   through pipeline. This will create "pixeldq" and "groupdq" arrays (if they
   don't already exist).

 - Propagate the DQ flags from the reference file DQ array to the science data "PIXELDQ"
   array using numpy's ``bitwise_or`` function.

Note that when applying the ``dq_init`` step to guide star data, the flags from the MASK reference
file are propagated into the guide star dataset "dq" array, instead of the "pixeldq" array.
The step identifies guide star data based on the following exposure type (exposure.type keyword attribute) values:
WFI_WIM_ACQ, WFI_WIM_TRACK, WFI_WSM_ACQ1, WFI_WSM_ACQ2, WFI_WSM_TRACK.

Conversion from Level 1 uncal files
-----------------------------------

The output of the initialization step is the RampModel. This is the form of the
data used throughout most of the exposure pipeline steps. For the most part, the
meta information between the input raw, or "uncal", model and the ramp model is
complete.

However, romancal supports processing a selection of files which use an
outdated schema. It supports these with a bespoke method that converts the files
to the new format when they are read in dq_init. This conversion does not do a
detailed mapping between all of the new and old metadata, but instead
opportunistically looks for fields with common names and assigns them. Other
metadata with non-matching names is simply copied in place. This allows
processing to proceed and preserves the original metadata, but the resulting
files have duplicates of many entries.
