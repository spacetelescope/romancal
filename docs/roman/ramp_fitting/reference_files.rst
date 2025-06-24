Reference Files
===============
The ``ramp_fit`` step uses three reference file types:

  - :ref:`GAIN <gain_reffile>`
  - :ref:`READNOISE <readnoise_reffile>`
  - :ref:`DARK <dark_reffile>`

During ramp fitting, the GAIN values are used to temporarily convert the pixel
values from units of DN to electrons, and convert the results of ramp fitting
back to DN.  The READNOISE values are used as part of the noise estimate for
each pixel. The DARK values are needed to estimate the Poisson noise. See
:ref:`rampfit-algorithm-ols` for more details.
