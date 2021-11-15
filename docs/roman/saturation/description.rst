Description
============

The ``saturation`` step flags pixels at or below the A/D floor or above the
saturation threshold.  Pixels values are flagged as saturated if the pixel value
is larger than the defined saturation threshold.  Pixel values are flagged as
below the A/D floor if they have a value of zero DN.

This step examines the data group-by-group, comparing the pixel values in the data array with defined
saturation thresholds for each pixel. When it finds a pixel value in a given
group that is above the saturation threshold (high saturation), it sets the
"SATURATED" flag in the corresponding location of the "groupdq" array in the
science exposure.  When it finds a pixel in a given group that has a zero or
negative value (below  the A/D floor), it sets the "AD_FLOOR" and "DO_NOT_USE"
flags in the corresponding location of the "groupdq" array in the science
exposure  For the saturation case, it also flags all subsequent groups for that
pixel as saturated. For example, if there are 10 groups and
group 7 is the first one to cross the saturation threshold for a given pixel,
then groups 7 through 10 will all be flagged for that pixel.

Pixels with thresholds set to NaN or flagged as "NO_SAT_CHECK" in the saturation
reference file have their thresholds set to the 16-bit A-to-D converter limit
of 65535 and hence will only be flagged as saturated if the pixel reaches that
hard limit in at least one group. The "NO_SAT_CHECK" flag is propagated to the
PIXELDQ array in the output science data to indicate which pixels fall into
this category.
