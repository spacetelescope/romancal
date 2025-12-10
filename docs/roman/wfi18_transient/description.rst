Description
============

The ``wfi18_transient`` step corrects for an anomalous transient signal in the first
read of exposures taken with detector WFI18.  The transient signal is strongly temperature
dependent with significant exposure-to-exposure variations so it must be fit and
removed separately for each exposure.

To fit the transient signal, the ``wfi_transient`` step implements the following algorithm:

#. Estimate the true counts in the first read from the values in the 2nd and 3rd
   resultants, extrapolating backward from the count rate.

#. Construct the first read residuals by subtracting the estimated counts from the
   measured counts, recorded in the first resultant.

#. Identify weakly illuminated pixels in each row, using the 4th and 5th resultant
   values to minimize bias from covariances.

#. Estimate the typical residual in each row from a sigma-clipped mean of the
   weakly illuminated pixels.

#. Fit the row residuals with a double-exponential model with functional form:

   .. math::

      residual = a\ exp[-t / \tau_a] + b\ exp[-t / \tau_b]

   where :math:`t` is the average readout time for the row and :math:`a`, :math:`b`,
   :math:`\tau_a` and :math:`\tau_b` are free parameters.

#. Evaluate the fit model at the readout time for each pixel in the first
   resultant.

#. Subtract the modeled residuals from the first resultant.

This algorithm requires a minimum of 5 resultants to be successful.  If fewer resultants
are present, or if the fit fails for any reason, the step will fall back on simply
masking the rows most affected by the transient anomaly.  In this case, the first
1000 rows of the first resultant will be set to DO_NOT_USE in the GROUPDQ array.
Optionally, the user may choose to apply the mask instead of attempting a fit as above,
by setting the ``mask_rows`` parameter for this step to ``True``.

This step is applicable to detector WFI18 only.  On successful completion for WFI18
data, the ``wfi18_transient`` keyword in the ``cal_step`` metadata is set to "COMPLETE".
If it is run on an exposure from any other detector it will have no effect, other
than to set the ``cal_step`` metadata keyword to "N/A".

References
----------

The WFI18 transient correction algorithm is based on work by T. Brandt,
"Characterizing and Correcting the Anomalous First Read of Detector 18 on Roman-WFI",
STScI Technical Document, 2025 (in prep).
