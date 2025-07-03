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

The current implementation uses fitted dark slope information from the dark reference file to
remove the dark current contribution from the WFI data.
The read pattern and frame time is retrieved from the input data file and the dark
reference file supplies the slope fitted to the dark current. This information is used to
predict the contribution of the dark current, or, more accurately, the background, as
it includes thermal emission from the telescope and instrument is subtracted from the
2-d image data produced by the ramp fitting step. 


The ERR arrays of the science data are updated and the error in the dark current
is added in quadrature to the existing error.  

The DQ flags from the dark reference file are propagated into science
exposure PIXELDQ array using a bitwise OR operation.

Upon successful completion of the dark subtraction the cal_step attribute is
set to "COMPLETE".

Special Handling
++++++++++++++++

Any pixel values in the dark reference data that are set to NaN will have their
values reset to zero before being subtracted from the science data, which
will effectively skip the dark subtraction operation for those pixels.
