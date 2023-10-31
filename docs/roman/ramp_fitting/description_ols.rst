.. _rampfit-algorithm-ols:

====================================================
Description: Optimized Least-squares with Even Ramps
====================================================

This step determines the mean count rate, in units of counts per second, for
each pixel by performing a linear fit to the data in the input file.  The fit
is done using the "ordinary least squares" method.
The fit is performed independently for each pixel.  There can be up to two
output files created by the step. The primary output file ("rate") contains the
slope at each pixel.
A second, optional output product is also available, containing detailed fit
information for each pixel. The two types of output files are described in
more detail below.

The count rate for each pixel is determined by a linear fit to the
cosmic-ray-free and saturation-free ramp intervals for each pixel; hereafter
this interval will be referred to as a "segment." The fitting algorithm uses an
'optimal' weighting scheme, as described by Fixsen et al, PASP, 112, 1350.
Segments are determined using
the 3-D GROUPDQ array of the input data set, under the assumption that the jump
step will have already flagged CR's. Segments are terminated where
saturation flags are found. Pixels are processed simultaneously in blocks
using the array-based functionality of numpy.  The size of the block depends
on the image size and the number of groups.

The ramp fitting step is also where the :ref:`reference pixels <refpix>` are
trimmed, resulting in a smaller array being passed to the subsequent steps.

Multiprocessing
===============

This step has the option of running in multiprocessing mode. In that mode it will
split the input data cube into a number of row slices based on the number of available
cores on the host computer and the value of the max_cores input parameter. By
default the step runs on a single processor. At the other extreme if max_cores is
set to 'all', it will use all available cores (real and virtual). Testing has shown
a reduction in the elapsed time for the step proportional to the number of real
cores used. Using the virtual cores also reduces the elasped time but at a slightly
lower rate than the real cores.

Special Cases
+++++++++++++

If the input dataset has only a single group, the count rate
for all unsaturated pixels will be calculated as the
value of the science data in that group divided by the group time.  If the
input dataset has only two groups, the count rate for all
unsaturated pixels will be calculated using the differences
between the two valid groups of the science data.

For datasets having more than a single group, a ramp having
a segment with only a single group is processed differently depending on the
number and size of the other segments in the ramp. If a ramp has only one
segment and that segment contains a single group, the count rate will be calculated
to be the value of the science data in that group divided by the group time.  If a ramp
has a segment having a single group, and at least one other segment having more
than one good group, only data from the segment(s) having more than a single
good group will be used to calculate the count rate.

The data are checked for ramps in which there is good data in the first group,
but all first differences for the ramp are undefined because the remainder of
the groups are either saturated or affected by cosmic rays.  For such ramps,
the first differences will be set to equal the data in the first group.  The
first difference is used to estimate the slope of the ramp, as explained in the
'segment-specific computations' section below.

If any input dataset contains ramps saturated in their second group, the count
rates for those pixels will be calculated as the value
of the science data in the first group divided by the group time.

All Cases
+++++++++
For all input datasets, including the special cases described above, arrays for
the primary output (rate) product are computed as follows.

After computing the slopes for all segments for a given pixel, the final slope is
determined as a weighted average from all segments, and is
written as the primary output product.  In this output product, the
3-D GROUPDQ is collapsed into 2-D, merged
(using a bitwise OR) with the input 2-D PIXELDQ, and stored as a 2-D DQ array.

A second, optional output product is also available and is produced only when
the step parameter 'save_opt' is True (the default is False).  This optional
product contains 3-D arrays called SLOPE, SIGSLOPE, YINT, SIGYINT, WEIGHTS,
VAR_POISSON, and VAR_RNOISE that contain the slopes, uncertainties in the
slopes, y-intercept, uncertainty in the y-intercept, fitting weights, the
variance of the slope due to poisson noise only, and the variance of the slope
due to read noise only for each segment of each pixel, respectively. The y-intercept refers
to the result of the fit at an effective exposure time of zero.  This product also
contains a 2-D array called PEDESTAL, which gives the signal at zero exposure
time for each pixel, and the 3-D CRMAG array, which contains the magnitude of
each group that was flagged as having a CR hit.  By default, the name of this
output file will have the suffix "_fitopt".
In this optional output product, the pedestal array is
calculated by extrapolating the final slope (the weighted
average of the slopes of all ramp segments) for each pixel
from its value at the first group to an exposure time of zero. Any pixel that is
saturated on the first group is given a pedestal value of 0. Before compression,
the cosmic ray magnitude array is equivalent to the input SCI array but with the
only nonzero values being those whose pixel locations are flagged in the input
GROUPDQ as cosmic ray hits. The array is compressed, removing all groups in
which all the values are 0 for pixels having at least one group with a non-zero
magnitude. The order of the cosmic rays within the ramp is preserved.

Slope and Variance Calculations
+++++++++++++++++++++++++++++++
Slopes and their variances are calculated for each segment,
and for the entire exposure. As defined above, a segment is a set of contiguous
groups where none of the groups are saturated or cosmic ray-affected.  The
appropriate slopes and variances are output to the primary output product, and the optional output product. The
following is a description of these computations. The notation in the equations
is the following: the type of noise (when appropriate) will appear as the superscript
‘R’, ‘P’, or ‘C’ for readnoise, Poisson noise, or combined, respectively;
and the form of the data will appear as the subscript: ‘s’, ‘o’ for segment, or overall (for the entire dataset), respectively.

Optimal Weighting Algorithm
---------------------------
The slope of each segment is calculated using the least-squares method with optimal
weighting, as described by Fixsen et al. 2000, PASP, 112, 1350; Regan 2007,
JWST-STScI-001212. Optimal weighting determines the relative weighting of each sample
when calculating the least-squares fit to the ramp. When the data have low signal-to-noise
ratio :math:`S`, the data are read noise dominated and equal weighting of samples is the
best approach. In the high signal-to-noise regime, data are Poisson-noise dominated and
the least-squares fit is calculated with the first and last samples. In most practical
cases, the data will fall somewhere in between, where the weighting is scaled between the
two extremes.

The signal-to-noise ratio :math:`S` used for weighting selection is calculated from the
last sample as:

.. math::
    S = \frac{data \times gain} { \sqrt{(read\_noise)^2 + (data \times gain) } } \,,

The weighting for a sample :math:`i` is given as:

.. math::
    w_i = (i - i_{midpoint})^P \,,

where :math:`i_{midpoint}` is the the sample number of the midpoint of the sequence, and
:math:`P` is the exponent applied to weights, determined by the value of :math:`S`. Fixsen
et al. 2000 found that defining a small number of P values to apply to values of S was
sufficient; they are given as:

+-------------------+------------------------+----------+
| Minimum S         | Maximum S              | P        |
+===================+========================+==========+
| 0                 | 5                      | 0        |
+-------------------+------------------------+----------+
| 5                 | 10                     | 0.4      |
+-------------------+------------------------+----------+
| 10                | 20                     | 1        |
+-------------------+------------------------+----------+
| 20                | 50                     | 3        |
+-------------------+------------------------+----------+
| 50                | 100                    | 6        |
+-------------------+------------------------+----------+
| 100               |                        | 10       |
+-------------------+------------------------+----------+

Segment-specific Computations:
------------------------------
The variance of the slope of a segment due to read noise is:

.. math::
   var^R_{s} = \frac{12 \ R^2 }{ (ngroups_{s}^3 - ngroups_{s})(group_time^2) } \,,

where :math:`R` is the noise in the difference between 2 frames,
:math:`ngroups_{s}` is the number of groups in the segment, and :math:`group_time` is the group
time in seconds (from the exposure.group_time).

The variance of the slope in a segment due to Poisson noise is:

.. math::
   var^P_{s} = \frac{ slope_{est} }{  tgroup \times gain\ (ngroups_{s} -1)}  \,,

where :math:`gain` is the gain for the pixel (from the GAIN reference file),
in e/DN. The :math:`slope_{est}` is an overall estimated slope of the pixel,
calculated by taking the median of the first differences of the groups that are
unaffected by saturation and cosmic rays. This is a more
robust estimate of the slope than the segment-specific slope, which may be noisy
for short segments.

The combined variance of the slope of a segment is the sum of the variances:

.. math::
   var^C_{s} = var^R_{s} + var^P_{s}


Exposure-level computations:
----------------------------

The variance of the slope due to read noise is:

.. math::
   var^R_{o} = \frac{1}{ \sum_{s} \frac{1}{ var^R_{s}}}

where the sum is over all segments.


The variance of the slope due to Poisson noise is:

.. math::
   var^P_{o} = \frac{1}{ \sum_{s} \frac{1}{ var^P_{s}}}

The combined variance of the slope is the sum of the variances:

.. math::
   var^C_{o} = var^R_{o} + var^P_{o}

The square root of the combined variance is stored in the ERR array of the primary output.

The overall slope depends on the slope and the combined variance of the slope of all
segments, so is a sum over segments:

.. math::
    slope_{o} = \frac{ \sum_{s}{ \frac{slope_{s}} {var^C_{s}}}} { \sum_{s}{ \frac{1} {var^C_{s}}}}


Upon successful completion of this step, the status attribute ramp_fit will be set
to "COMPLETE".


Error Propagation
=================

Error propagation in the ramp fitting step is implemented by storing the
square-root of the exposure-level combined variance in the ERR array of the primary
output product. This combined variance of the exposure-level slope is the sum
of the variance of the slope due to the Poisson noise and the variance of the
slope due to the read noise. These two variances are also separately written
to the arrays VAR_POISSON and VAR_RNOISE in the asdf output.

For the optional output product, the variance of the slope due to the Poisson
noise of the segment-specific slope is written to the VAR_POISSON array.
Similarly, the variance of the slope due to the read noise of the
segment-specific slope  is written to the VAR_RNOISE array.
