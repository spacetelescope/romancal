Description
============

Overview
--------

The ``refpix`` step corrects for additional signal from the electronics using
the reference pixels.

Reference Pixels in Data
------------------------

The WFI has two sets of reference pixels: a 4-pixel border of reference pixels
around the science pixels, and the Amp33 reference pixels which are a
4096 x 128 section of the detector adjacent to the border pixels.

In the data files, the storage location of the reference pixels depends on level
of processing.

A Level 1, uncalibrated image has one array that contains both the science
pixels and the border reference pixel, and a separate array for the Amp33 pixels.

A RampModel, which is created during the ``dq_init`` step, represents a dataset
at any intermediate step between Level 1 and the final Level 2 image. Like the
Level 1 file, RampModels also contain an array with both the science and border
reference pixels, and another with the Amp 33 reference pixels. In addition to
these arrays, there are four more arrays that contain the original border
reference pixels (top, bottom, left, and right), and an additional four for
their DQ arrays. The border pixels are copied during the ``dq_init``, so they
reflect the original state of the border pixels before any calibration.
The border pixels that are still attached to the science data in the RampModel
will later be discarded when the Level 2 image is created. Note that the border
reference pixel arrays each include the overlap regions in the corners, so that
each slice contains the full span of border pixels at the top, bottom, left, or
right.

In the Level 2, calibrated image, the data array only contains the science
pixels. The border reference pixels are trimmed from this image duing
ramp_fit. The additional arrays for the original border reference pixels
(which are 3D) and their DQ arrays, and the Amp 33 reference pixels, remain in
the Level 2 file.

Discretization bias & reference pixel correction
------------------------------------------------

The analog-to-digital conversion in the Roman electronics performs an
integer floor operation that biases the downlinked signal low relative
to the actual number of photons observed by the instrument.  The
equation for this "discretization bias" is given by:

.. math:: \mathrm{bias} = -0.5 - 0.5 \frac{N-1}{N} \, ,

in units of counts, where :math:`N` is the number of reads entering
into a particular resultant.  This is a small effect.  The constant :math:`-0.5`
term is degenerate with the pedestal and has no effect on ramp slopes
and therefore on the primary astronomical quantity of interest.  The
second term, however, depends on the number of reads in a resultant
and may vary from resultant to resultant in Roman.  This, if
uncorrected, can lead to a bias in the fluxes we derive from Roman
data for sources.

However, we need take no special action to correct for this effect.
The reference pixels are affected by the discretization bias in the
same way as the science pixels, and so when the reference pixels are
subtracted (roughly speaking!) from the science pixels, this bias cancels.
Exactly when this cancellation occurs depends on the details of the reference
pixel correction step.  Presently the reference pixel correct includes
a component that removes trends across each amplifier and frame using
the reference pixels at the top and bottom of the amplifier.  This
removes the discretization bias.

We note that even if the discretization bias were not removed at the
reference pixel correction stage, it could be corrected at the dark
subtraction step.  Provided that dark reference images are processed
through the usual reference pixel correction step, they will have the
same biases present in the reference-pixel-corrected images.  We have
decided to perform the dark subtraction of Roman images via
subtracting precomputed images for each MA table rather than scaling a
fixed dark rate image by the mean time of each resultant.  These
precomputed dark images will contain not only the dark current but
also electronic effects like the discretization bias.
However, it is better to correct this effect during the
reference pixel correction so that the dark reference images better
represent the dark current and can be more easily used to compute
Poisson uncertainties stemming from dark current.
