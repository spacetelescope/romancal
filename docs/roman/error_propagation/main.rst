Description
-----------
Steps in the various pipeline modules calculate variances due to different sources of
noise or modify variances that were computed by previous steps.
For some cases these arrays are being propagated to subsequent steps in the pipeline.
Anytime a step creates or updates variances, the total error (ERR) array values
are recomputed as the square root of the quadratic sum of all variances available
to the step.
Note that the ERR array values are always expressed as a standard deviation
(the square root of the variance).

The table below is a summary of which steps create or update variance and error arrays,
as well as which steps make use of these data. Details of how each step computes or
uses these data are given in the subsequent sections below.

================= ===== ======================= ============================ =========
Step              Stage Creates arrays          Updates arrays               Step uses
================= ===== ======================= ============================ =========
ramp_fitting      ELPP  VAR_POISSON, VAR_RNOISE ERR                          None
flat_field        ELPP  VAR_FLAT                ERR, VAR_POISSON, VAR_RNOISE None
outlier_detection HLPP  None                    None                         ERR
================= ===== ======================= ============================ =========

ELPP Processing
---------------
ELPP processing pipelines perform detector-level corrections and ramp fitting for
individual exposures, for nearly all imaging and spectroscopic modes. Details
of the pipelines can be found at :ref:`roman_elp <exposure_pipeline>`.

The ELPP pipeline steps that alter the ERR, VAR_POISSON, or VAR_RNOISE arrays
the science countrate data are discussed below.
Any step not listed here does not alter or use the variance or error arrays
in any way and simply propagates the information to the next step.

ramp_fitting
++++++++++++
This step calculates and populates the VAR_POISSON and VAR_RNOISE arrays to pass to the
next step or saved in the optional output 'rate' files. The ERR array is updated as the square root of the
quadratic sum of the variances. VAR_POISSON and VAR_RNOISE represent the uncertainty in the
computed slopes (per pixel) due to Poisson and read noise, respectively.
The details of the calculations can be found at :ref:`ramp_fitting <ramp_fitting_step>`.

flat_field
++++++++++
The SCI array of the exposure being processed is divided by the flat-field reference
image, and the VAR_POISSON and VAR_RNOISE arrays are divided by the square of the flat.
A VAR_FLAT array is created, computed from the science data and the flat-field
reference file ERR array.
The science data ERR array is then updated to be the square root of the quadratic sum of
the three variance arrays.
For more details see :ref:`flat_field <flatfield_step>`.

HLPP Processing
---------------
HLPP pipelines perform additional instrument-level and observing-mode corrections and
calibrations to produce fully calibrated exposures. The various steps that apply corrections and calibrations apply those same corrections/calibrations to all variance arrays and update the total
ERR array.

outlier_detection
+++++++++++++++++
The ``outlier_detection`` step is used in all Stage 3 pipelines.
