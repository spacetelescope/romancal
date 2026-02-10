Description
============

Assumptions
-----------
It is assumed that the saturation step has already been applied to
the input data, so that saturation flags are set in the GROUPDQ array of
the input science data.

Algorithm
---------
The linearity step corrects science data values for detector non-linearity.
The correction is applied pixel-by-pixel, group-by-group within a science
exposure, using the linearity correction routines in stcal.

The correction has two components:

1. **Integral Nonlinearity (INL)**: A per-channel lookup table correction that
   accounts for integral nonlinearity in the analog-to-digital converter.

2. **Classical Nonlinearity (CNL)**: A polynomial correction that accounts for
   the change in pixel capacitance with accumulated charge.

Since the pipeline operates on resultants (averages of individual reads) rather
than individual reads, the stcal routines simulate the individual reads from
the observed resultants based on the read pattern. This simulation requires
inverse linearity polynomials. The corrections are applied to each simulated
read, and then the reads are re-averaged into corrected resultants.

Upon successful completion of the linearity correction, "cal_step" in the
metadata is set to "COMPLETE".

Classical Nonlinearity Correction
+++++++++++++++++++++++++++++++++
The classical nonlinearity correction is represented by an nth-order polynomial
for each pixel in the detector, with n+1 arrays of coefficients read from the
linearity reference file.

The algorithm for correcting the observed pixel value in each group is:

.. math::
   F_\text{c} = c_{0} + c_{1}F + c_{2}F^2 + c_{3}F^3 + ... + c_{n}F^n

where :math:`F` is the observed counts (in DN), :math:`c_n` are the polynomial
coefficients, and :math:`F_\text{c}` is the corrected counts. There is no
limit to the order of the polynomial correction; all coefficients contained in
the reference file will be applied.

Integral Nonlinearity Correction
++++++++++++++++++++++++++++++++
The integral nonlinearity (INL) correction addresses nonlinearity in the
analog-to-digital converter. Each of the 32 readout channels (128 columns each)
has its own correction table that maps observed DN values to a correction
offset. See `Brandt & Perera (2025)
<https://www.stsci.edu/files/live/sites/www/files/home/roman/documentation/technical-documentation/_documents/Roman-STScI–000866.pdf>`_
for details.

The INL correction is optional; if no integral nonlinearity reference file is
available, only the classical polynomial correction is applied.

Special Handling
++++++++++++++++

- Pixels having at least one correction coefficient equal to NaN will not have
  the linearity correction applied and the DQ flag "NO_LIN_CORR" is added to
  the science exposure PIXELDQ array.

- Pixels that have the "NO_LIN_CORR" flag set in the DQ array of the linearity
  reference file will not have the correction applied and the "NO_LIN_CORR" flag
  is added to the science exposure PIXELDQ array.

- Pixel values that have the "SATURATED" flag set in a particular group of the
  science exposure GROUPDQ array will not have the linearity correction
  applied to that group. Any groups for that pixel that are not flagged as
  saturated will be corrected.

The ERR array of the input science exposure is not modified.

The flags from the linearity reference file DQ array are propagated into the
PIXELDQ array of the science exposure using a bitwise OR operation.
