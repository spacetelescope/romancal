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
dark current rate data stored in a dark reference file from the fitted ramp
data. 

The current implementation uses dark reference files for a given
detector that have the rate of the accumulated dark current measured
and stored in the `dark_slope` array in the reference file.
The dark rate are then subtracted from the fitted ramp data. 

The reference dark rate will have an associated uncertainty and this error is added
in quadrature to the error determined in the ramp fit step. 

The DQ flags from the dark reference file are propagated into science
exposure DQ array using a bitwise OR operation.

Upon successful completion of the dark subtraction the cal_step attribute is
set to "COMPLETE".

Special Handling
++++++++++++++++

Any pixel values in the dark reference data that are set to NaN will have their
values reset to zero before being subtracted from the science data, which
will effectively skip the dark subtraction operation for those pixels.
