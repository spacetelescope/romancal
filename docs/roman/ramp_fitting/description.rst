.. _rampfit-algorithm-ols22:

Description: Optimized Least-squares with Uneven Ramps
======================================================

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
this interval will be referred to as a "segment." There are two algorithms used:
Optimal Least-Square ('ols') and Optimal Least-Square for Uneven ramps
('ols_cas22'). The 'ols' algorithm is the one
`used by JWST <https://jwst-pipeline.readthedocs.io/en/stable/jwst/ramp_fitting/description.html>`_
and is further :ref:`described here <rampfit-algorithm-ols>`.

The 'ols_22' algorithm is based on `Casertano et al, STScI Technical Document,
2022
<https://www.stsci.edu/files/live/sites/www/files/home/roman/_documents/Roman-STScI-000394_DeterminingTheBestFittingSlope.pdf>`_.
The implementation is what is described in this document.

Segments
++++++++

Segments are determined using the 3-D GROUPDQ array of the input data set, under
the assumption that the jump step will have already flagged CR's. Segments are
terminated where saturation flags are found.

The ramp fitting step is also where the :ref:`reference pixels <refpix>` are
trimmed, resulting in a smaller array being passed to the subsequent steps.

Special Cases
+++++++++++++

If the input dataset has only a single resultant, no fit is determined, giving
that resultant a weight of zero.

All Cases
+++++++++
For all input datasets, including the special cases described above, arrays for
the primary output (rate) product are computed as follows.

After computing the slopes for all segments for a given pixel, the final slope is
determined as a weighted average from all segments, and is
written as the primary output product.  In this output product, the
3-D GROUPDQ is collapsed into 2-D, merged
(using a bitwise OR) with the input 2-D PIXELDQ, and stored as a 2-D DQ array.

Slope and Variance Calculations
+++++++++++++++++++++++++++++++
Slopes and their variances are calculated for each segment,
and for the entire exposure. As defined above, a segment is a set of contiguous
resultants where none of the resultants are saturated or cosmic ray-affected.  The
appropriate slopes and variances are output to the primary output product, and the optional output product. The
following is a description of these computations. The notation in the equations
is the following: the type of noise (when appropriate) will appear as the superscript
‘R’, ‘P’, or ‘C’ for readnoise, Poisson noise, or combined, respectively;
and the form of the data will appear as the subscript: ‘s’, ‘o’ for segment, or overall (for the entire dataset), respectively.

Optimal Weighting Algorithm
---------------------------
The slope of each segment is calculated using the least-squares method with optimal
weighting, as described by `Casertano et al, STScI Technical Document,
2022
<https://www.stsci.edu/files/live/sites/www/files/home/roman/_documents/Roman-STScI-000394_DeterminingTheBestFittingSlope.pdf>`_.
Optimal weighting determines the relative weighting of each sample
when calculating the least-squares fit to the ramp. When the data have low signal-to-noise
ratio :math:`S`, the data are read noise dominated and equal weighting of samples is the
best approach. In the high signal-to-noise regime, data are Poisson-noise dominated and
the least-squares fit is calculated with the first and last samples. In most practical
cases, the data will fall somewhere in between, where the weighting is scaled between the
two extremes.

The signal-to-noise ratio :math:`S` used for weighting selection is calculated from the
last sample as:

.. math::
   S_{max} = S_{last} - S_{first}

   S = \frac{S_{max}} { \sqrt{(read\_noise)^2 + S_{max} } } \,,

where :math:`S_{max}` is the maximum signal in electrons with the pedestal
removed.

The weighting for a sample :math:`i` is given as:

.. math::
    w_i = \frac{(1 + P) \times N_i} {1 + P \times N_i} | \bar t_i - \bar t_{mid} |^P \,,

where :math:`t_{mid}` is the time midpoint of the sequence,
:math:`N_i` is the number of contributing reads, and
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

Segment-Specific Computations
-----------------------------

The segment fitting implementation is based on Section 5 of Casertano et.
al. 2022. Segments which have only a single resultant, no fitting is performed.

A set of auxiliary quantities are computed as follows:

.. math::
   F0 &= \sum_{i=0}^{N-1} W_i

   F1 &= \sum_{i=0}^{N-1} W_i \bar t_i

   F2 &= \sum_{i=0}^{N-1} W_i \bar t_i^2

The denominator, :math:`D`, is calculated as a single two-dimensional array:

.. math::
   D = F2 \cdot F0 - F1^2


The resultant coefficients, :math:`K_i`, are computed as N two dimensional
arrays:

.. math::
   K_i = (F0 \cdot \bar t_i - F1) \cdot W_i / D

The estimated slope, :math:`\hat F`, is computed as a sum over the resultants
:math:`R_i` and the coefficients :math:`K_i`:

.. math::
   \hat F = \sum_{i} K_i R_i

The read-noise component :math:`V_R` of the slope variance is computed as:

.. math::
   V_R = \sum_{i=0}^{N-1} K_i^2 \cdot (RN)^2 / N_i

The signal variance, :math:`V_S`, of the count rate in the signal-based component of the slope
variance is computed as:

.. math::
   V_S = \sum_{i=0}^{N-1} {K_i^2 \tau_i} + \sum_{i<j} {2 K_i K_j \cdot \bar t_i}

Total variance, if desired, is a estimate of the total slope variance :math:`V` can
be computed by adopting :math:`\hat F` as the estimate of the slope:

.. math::
   V = V_R + V_S \cdot \hat F

Exposure-level computations:
----------------------------

The ramps for each resultant are reconstructed from its segments, :math:`i`,
fits by calculating the inverse variance-weighted mean using the read noise
variances:

.. math::
   w_i &= 1 / V_{R_i}

   \hat F_{mean} &= \frac {\sum_i {w_i \hat F_i}} {\sum_i w_i}

The read noise is determined as follows:

.. math::
   V_{R_{mean}} = \frac {\sum_i {w_i ^ 2 V_{R_i}}} {(\sum_i {w_i}) ^ 2}

Finally, the signal variance is calculated as:

.. math::

   V_{S_{mean}} = \frac {\sum_i {w_i ^ 2 V_{S_i}}} {(\sum_i {w_i}) ^ 2}

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
