Description
============

The ``gain_scale`` step rescales pixel values in roman countrate
science data products in order to correct for the effect of using
a non-standard detector gain setting. The countrate data are
rescaled to make them appear as if they had been obtained using
the standard gain setting.

The `gain_scale` step is applied at the end of the
TBD 
#:ref:`calwebb_detector1 <calwebb_detector1>` 
pipeline, after the ramp_fit <ramp_fitting_step> step has been applied. It is applied to both the "rate" and "rateints" products from
ramp_fit <ramp_fitting_step>, if both
types of products were created. The science (SCI) and error (ERR)
arrays are multiplied by the gain factor, and the Poisson
variance (VAR_POISSON) and read noise variance (VAR_RNOISE) arrays
are multiplied by the square of the gain factor.

The scaling factor is obtained from the "gain_factor" keyword attribute in the
header of the gain reference file. Normally the
ramp_fit <ramp_fitting_step> step
reads that keyword attribute value during its execution and stores the value in
the science data "gain_factor" keyword attribute, so that the gain reference file
does not have to be loaded again by the ``gain_scale`` step. If, however,
the step does not find that keyword attribute populated in the science data, it
loads the gain reference file to retrieve it. If all attempts to
find the scaling factor fail, the step is skipped.

Gain reference files for instruments or modes that use the standard
gain setting will typically not have the "gain_factor" keyword attribute in their
metadata, which causes the ``gain_scale`` step to be skipped. Alternatively,
gain reference files for modes that use the standard gain can have
gain_factor=1.0, in which case the correction is benign.

Upon successful completion of the step, the <TBD> keyword attribute in the
science data is set to "COMPLETE".
