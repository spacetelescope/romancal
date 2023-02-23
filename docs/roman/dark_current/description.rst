Description
===========

Assumptions
-----------

It is assumed that the input science data have had the zero group (or
bias) subtracted. Accordingly, the dark reference data
should have their own group zero subtracted from all groups.

Algorithm
---------

The dark current step removes dark current from a Roman exposure by subtracting
dark current data stored in a dark reference file.

The current implementation uses dark reference files that are matched to the
MA table entry in the exposure metadata. Note that the data reference file
for a science group (SCI) is named `data`. The dark data are then subtracted,
group-by-group, from the science exposure groups, in which
each SCI group of the dark data is subtracted from the corresponding SCI
group of the science data.

The ERR arrays of the science data are not modified.

The DQ flags from the dark reference file are propagated into science
exposure PIXELDQ array using a bitwise OR operation.

Upon successful completion of the dark subtraction the cal_step attribute is
set to "COMPLETE".

Special Handling
++++++++++++++++

Any pixel values in the dark reference data that are set to NaN will have their
values reset to zero before being subtracted from the science data, which
will effectively skip the dark subtraction operation for those pixels.
