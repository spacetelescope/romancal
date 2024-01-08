.. _highlevel_pipeline:


High Level Image Processing
=====================================================

:Class: `romancal.pipeline.HighLevelPipeline`
:Alias: highlevel_pipeline

The ``HighLevelPipeline`` applies corrections to an overlapping group of images
and is setup to process only imaging observations.
This pipeline is used to determine a common background, :ref:`skymatch <skymatch_step>`, detect pixels the are
not consistent with the other datasets, :ref:`outlier_detection <outlier_detection_step>`, and resample the image to a
single undistorted image, :ref:`resample <resample_step>`.

The list of steps applied by the ``HighLevelPipeline`` pipeline is shown in the
table below.

.. |check| unicode:: U+2713 .. checkmark
.. |xmark| unicode:: U+1D54F .. xmark

======================================================= ========= ========= =========
 Step                                                   WFI-Image WFI-Prism WFI-Grism
======================================================= ========= ========= =========
 :ref:`skymatch <skymatch_step>`                        |check|    |xmark|  |xmark|
 :ref:`outlier_detection <outlier_detection_step>`      |check|    |xmark|  |xmark|
 :ref:`resample <resample_step>`                        |check|    |xmark|  |xmark|
======================================================= ========= ========= =========


Arguments
---------
The ``highlevel`` pipeline has no optional arguments:


You can see the options for strun using:

strun --help roman_hlp

and this will list all the strun options all well as the step options for the roman_hlp.


Inputs
--------

2D image data
+++++++++++++

:Data model: `~romancal.datamodels.WfiImage`
:File suffix: _cal

The input to the ``HighLevelPipeline`` is a group of calibrated exposures,
e.g. "r0008308002010007027_06311_0019_WFI01_cal.asdf", which contains the
calibrated data for the the exposures. The most convenient way to pass the list of
exposures to be processed with the high level pipeline is to use an association.
Instructions on how to create an input association an be found at :ref:`asn-from-list`.


Outputs
----------

2D Image model
++++++++++++++

:Data model: `~romancal.datamodels.WfiMosaic`
:File suffix: _i2d

Result of applying all the high level pipeline steps up through the
:ref:`resample <resample_step>` step is to produce data background corrected
and cleaned of outliers and resampled to a distortion free grid.
This is 2D image data, with additional attributes for the mosaicing information.
