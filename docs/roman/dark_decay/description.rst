Description
===========


The ``dark_decay`` step corrects for an anomalous, decaying signal over the first
few resultants of WFI exposures. This has a typical amplitude of 0.4 DN
and a time constant of 23 s.  The effect is modeled as a simple exponential
decay, with a time constant and amplitude that varies from detector to detector, and
which affects all pixels on each detector in the same way.  As the detector reads out,
time passes and the contribution of the dark decay signal becomes steadily smaller.

This step works by using the ``darkdecaysignal`` reference files in CRDS to obtain the
measured amplitude and time constant for a given detector.  It loops over resultants,
and determines when each pixel in each frame in a resultant was read out.  It evaluates
the dark decay signal at that time.  It then averages the dark decay signal over the
frames entering a resultant, producing a correction image for each resultant.  It
subtracts this correction from the measured pixels.

References
----------

The dark decay correction is based on work by S. Betti et al,  See `The Statistical Properties of Dark Ramps for
the Roman-WFI Detectors (STScI Technical Document, 2025) <https://www.stsci.edu/files/live/sites/www/files/home/roman/documentation/technical-documentation/_documents/Roman-STScI-000814.pdf>`_
for more information.  
