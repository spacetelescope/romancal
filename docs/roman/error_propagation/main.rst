Description
-----------
Steps in the various pipeline modules calculate variances due to different sources of
noise or modify variances that were computed by previous steps.  In some cases the
variance arrays are only used internally within a given step.  For several steps,
these arrays must be propagated to subsequent steps in the pipeline. Any time a step
creates or updates variances, the total error (ERR) array values are always recomputed
as the square root of the quadratic sum of all variances available at the time.
Note that the ERR array values are always expressed as standard deviation
(i.e. square root of the variance).

The table below is a summary of which steps create or update variance and error arrays,
as well as which steps make use of these data. Details of how each step computes or
uses these data are given in the subsequent sections below.

================= ===== ======================= ====================================== =========
Step              Stage Creates arrays          Updates arrays                         Step uses
================= ===== ======================= ====================================== =========
ramp_fitting        1   VAR_POISSON, VAR_RNOISE ERR                                    None
gain_scale          1   None                    ERR, VAR_POISSON, VAR_RNOISE           None
flat_field          2   VAR_FLAT                ERR, VAR_POISSON, VAR_RNOISE           None
photom              2   None                    ERR, VAR_POISSON, VAR_RNOISE, VAR_FLAT None
================= ===== ======================= ====================================== =========

Stage 1 Pipelines 
-----------------
Stage 1 pipelines perform detector-level corrections and ramp fitting for
individual exposures, for nearly all imaging and spectroscopic modes. Details 
of the pipelines can be found at :ref:`Stage 1 Pipelines <calroman_detector1>`.

The Stage 1 pipeline steps that alter the ERR, VAR_POISSON, or VAR_RNOISE arrays of
the science countrate data are discussed below.
Any step not listed here does not alter or use the variance or error arrays
in any way and simply propagates the information to the next step.

ramp_fitting
++++++++++++
This step calculates and populates the VAR_POISSON and VAR_RNOISE arrays
in the 'rate' and 'rateints' files, and updates the ERR array as the square root of the
quadratic sum of the variances. VAR_POISSON and VAR_RNOISE represent the uncertainty in the
computed slopes (per pixel) due to Poisson and read noise, respectively.
The details of the calculations can be found at :ref:`ramp_fitting <ramp_fitting_step>`.

gain_scale
++++++++++
The ``gain_scale`` step is applied after ``ramp_fitting``, and applies to both the 
rate and rateints products. The gain correction is applied to the ERR, 
VAR_POISSON, and VAR_RNOISE arrays.  The SCI and ERR arrays are multiplied by the
gain correction factor, and the variance arrays are multiplied by the square of
the gain correction factor. More details can be
found at :ref:`gain_scale <gain_scale_step>`.

Stage 2 Pipelines 
-----------------
Stage 2 pipelines perform additional instrument-level and observing-mode corrections and 
calibrations to produce fully calibrated exposures. There are two main Stage 2 pipelines:
one for imaging :ref:`calroman_image2 <calroman_image2>` and one for 
spectroscopy :ref:`calroman_spec2 <calroman_spec2>`.
In these pipelines, the various steps that apply corrections and calibrations
apply those same corrections/calibrations to all variance arrays and update the total
ERR array.

flat_field
++++++++++
The SCI array of the exposure being processed is divided by the flat-field reference
image, and the VAR_POISSON and VAR_RNOISE arrays are divided by the square of the flat.
A VAR_FLAT array is created, computed from the science data and the flat-field
reference file ERR array.
The science data ERR array is then updated to be the square root of the quadratic sum of
the three variance arrays.
For more details see :ref:`flat_field <flatfield_step>`.

photom
++++++ 
The calibration information for the ``photom`` step includes a scalar flux conversion
constant, as well as optional arrays of wavelength and relative response (as a
function of wavelength). The combination of the scalar conversion factor and any 2-D
response values is applied to the science data, including the SCI and ERR arrays,
as well as the variance (VAR_POISSON, VAR_RNOISE, and VAR_FLAT) arrays. The flux
calibration values are multiplied into the science exposure SCI and ERR arrays,
and the square of the calibration values is multiplied into all variance arrays.
For details of the photom correction, see :ref:`photom <photom_step>`.

