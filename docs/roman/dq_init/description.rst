Description
============
The Data Quality (DQ) initialization step in the calibration pipeline
populates the DQ mask for the input dataset. Flag values from the
appropriate static mask ("MASK") reference file in CRDS are copied into the
"PIXELDQ" array of the input dataset, because it is assumed that flags in the
mask reference file pertain to problem conditions that affect all groups and
integrations for a given pixel.

The actual process consists of the following steps:

 - Determine what MASK reference file to use via the interface to the bestref
   utility in CRDS.

 - Propagate the DQ flags from the reference file DQ array to the science data "PIXELDQ"
   array using numpy's ``bitwise_or`` function.

Note that when applying the ``dq_init`` step to guide star data, as is done in
the :ref:`calwebb_guider <calwebb_guider>` pipeline, the flags from the MASK reference
file are propagated into the guide star dataset "DQ" array, instead of the "PIXELDQ" array.
The step identifies guide star data based on the following exposure type (exposure.type keyword attribute) values:
WFI_WIM_ACQ, WFI_WIM_TRACK, WFI_WSM_ACQ1, WFI_WSM_ACQ2, WFI_WSM_TRACK.
