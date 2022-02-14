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