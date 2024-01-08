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


Jump Detection
==============

For most pixels, the ramp steadily accumulates flux from the sky as an integration
proceeds.  However, in relatively rare cases, a cosmic ray can pass through the
detector which instantaneously deposits a large amount of charge in a pixel.
This leads the resulting ramp to have a discontinuous *jump* in a particular read,
and accordingly to discontinuities in the resultants downlinked from the telescope.
The jump detection algorithm attempts to identify uncontaminated segments of ramps
for ramp fitting, so that the underlying astronomical signal can be extracted without
contamination from these jumps.

If the uneven-ramp jump detection algorithm is turned on (the default), the ramp fitting
algorithm is then run iteratively on a "queue" (list) of ramps. The queue is initialized
with the ramp(s).
Then following the algorithm presented in Sharma et al (2023) (in preparation),
the jump detection algorithm picks a ramp, say :math:`[R_0, \dots, R_M]`, out of the
queue and runs the ramp fitting algorithm on it. It then checks the resulting ramp for jumps.
If a jump is detected, then two sub-ramps are created from the passed in ramp which exclude
the resultants predicted to be affected by the jump. These sub-ramps are then added to the
queue. This process continues until the queue is empty.

.. note::

   It may not always be possible to create two sub-ramps around the resultants predicted to
   be part of a jump. For example if these jump resultants include the first, second, second-to-last,
   or last resultant of the ramp then it is not possible to create two meaningful sub-ramps, as one
   cannot run the ramp fitting algorithm on a ramp with zero or only one resultant. Therefore, in
   these cases, only one sub-ramp is created and added to the queue.

The method use for determining if and where a jump occurs is divided into two parts. First,
a *statistic*, :math:`S` and possible jump resultants are determined for the fitted ramp.
Then the statistic is compared against a threshold function, :math:`S_{\text{threshold}}` to determine
if a jump has occurred.


Statistic and Possible Jump Resultants
++++++++++++++++++++++++++++++++++++++

The statistic used to determine if a jump has occurred in the ramp, :math:`[R_0, \dots, R_M]`,
is computed from a list of statistics computed for each *single* and *double-difference* of
the resultants in the ramp.  By single-difference we mean the difference between two adjacent
resultants in the ramp, while double-difference refers to the difference between a resultant
and a resultant two steps away (the resultant adjacent to a resultant adjacent to the resultant
in question).

To compute these statistics, the single-difference excess slope :math:`\delta_{i, i+1}` and
the double-difference excess slope :math:`\delta_{i, i+2}` are computed as:

.. math::

   \delta_{i, i+1} &= \frac{R_{i+1} - R_i} {\bar t_{i+1} - \bar t_i} - \hat \alpha

   \delta_{i, i+2} &= \frac{R_{i+2} - R_i} {\bar t_{i+2} - \bar t_i} - \hat \alpha

where :math:`\hat \alpha` is the slope computed by the ramp fitting algorithm. The
variance in the excess slope:

.. math::

   Var(\delta_{i, j}) &= \frac {Var(R_j - R_i)} {(\bar t_j - \bar t_i)^2} + f_{corr}(\hat \alpha)

   Var(R_j - R_i) &= \sigma_{RN} \left( \frac{1}{N_j} + \frac{1}{N_i} \right) + \hat \alpha \left[\tau_j + \tau_i - \min(\bar t_j, \bar t_i) \right]

   f_{corr}(\hat \alpha) &= - \frac{\hat \alpha}{t_{M - 1} - t_0}

where :math:`\sigma_{RN}` is the read noise. The single-difference statistic, :math:`s_i^\prime`,
and double-difference statistic, :math:`s_i^{\prime\prime}` are,

.. math::

   s_i^\prime &= \frac{\delta_{i, i+1}} {\sqrt{Var(\delta_{i, i+1})}}

   s_i^{\prime\prime} &= \frac{\delta_{i, i+2}} {\sqrt{Var(\delta_{i, i+2})}}.

The statistic :math:`s_i` for each resultants :math:`0 \leq i \leq M - 1` (no differences from the last
resultant are possible) is then computed as:

.. math::
   :nowrap:

   \[
   s_i =
   \begin{cases}
   s_i^\prime & \text{if } i = M - 2\\
   \max(s_i^\prime, s_i^{\prime\prime}) & \text{otherwise}
   \end{cases}
   \]


Finally, :math:`S = \max(s_i)` is the statistic used to determine if a jump has occurred in the fitted
ramp. The possible jump resultants for this ramp are the resultants :math:`R_i` and :math:`R_{i+1}`,
where :math:`i = \arg\max(s_i)`.

Two possible jump resultants are necessary, because the statistics cannot determine which of the two
adjacent resultants is the one affected by the jump. This is because if the jump occurs near the last
read making up :math:`R_i`, then it might appear that :math:`R_{i+1}` has the jump, this jump will be
picked up the :math:`s_i^{\prime\prime}` statistic. Using just the :math:`s_i^\prime` statistic, the
jump would be incorrectly identified in :math:`R_{i+1}`.


Threshold Function
++++++++++++++++++

Similarly to the statistic, the threshold function relies on the slope computed by the ramp fitting
algorithm, :math:`\hat \alpha`. The function itself was determined empirically by running simulations
of ramps with jumps and ramps without jumps. The threshold function was determined to be:

.. math::

   S_{\text{threshold}}(\hat \alpha) = 5.5 - \frac{1}{3}\log_{10}(\hat \alpha)

This corresponds to identifying jumps at 5.5 sigma when the count rate is 1 electron per second, and
4.5 sigma when the count rate is 1000 electrons per second. The decision was made to have the threshold
depend on the count rate because the pixels with lots of signal have larger uncertainties; meaning that
lower amplitude cosmic rays get identified in these cases.


A jump is determined to have occurred for a ramp fit with statistic, :math:`S`, with possible jump
resultants :math:`R_i,\ R_{i+1}`, if :math:`S \geq S_{\text{threshold}}(\hat \alpha)`. This results
in the ramp being split into two sub-ramps :math:`[R_0, \dots R_{i-1}]` and :math:`[R_{i+2}, \dots R_M]`,
which are then added to the ramp queue.

Error Propagation
=================

Error propagation in the ramp fitting step is implemented by storing the
square-root of the exposure-level combined variance in the ERR array of the primary
output product. This combined variance of the exposure-level slope is the sum
of the variance of the slope due to the Poisson noise and the variance of the
slope due to the read noise. These two variances are also separately written
to the arrays VAR_POISSON and VAR_RNOISE in the asdf output.
