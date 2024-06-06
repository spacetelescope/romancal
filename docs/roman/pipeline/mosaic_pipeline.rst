.. _mosaic_pipeline:


Mosaic Level Image Processing
=============================

:Class: `romancal.pipeline.MosaicPipeline`
:Alias: mosaic_pipeline

The ``MosaicPipeline`` applies corrections to an overlapping group of images
and is setup to process only imaging observations.
This pipeline is used to determine a common background, :ref:`skymatch <skymatch_step>`, detect pixels the are
not consistent with the other datasets, :ref:`outlier_detection <outlier_detection_step>`, and resample the image to a
single undistorted image, :ref:`resample <resample_step>`.

The list of steps applied by the ``MosaicPipeline`` pipeline is shown in the
table below.

.. |check| unicode:: U+2713 .. checkmark
.. |xmark| unicode:: U+1D54F .. xmark

======================================================= ========= ========= =========
 Step                                                   WFI-Image WFI-Prism WFI-Grism
======================================================= ========= ========= =========
 :ref:`skymatch <skymatch_step>`                        |check|    |xmark|  |xmark|
 :ref:`outlier_detection <outlier_detection_step>`      |check|    |xmark|  |xmark|
 :ref:`resample <resample_step>`                        |check|    |xmark|  |xmark|
 :ref:`source_catalog <source_catalog_step>`            |check|    |xmark|  |xmark|
======================================================= ========= ========= =========


Arguments
---------
The ``mosaic`` pipeline has no optional arguments:


You can see the options for strun using:

strun --help roman_mos

and this will list all the strun options all well as the step options for the roman_mos.


Inputs
--------

An association of 2D calibrated image data
++++++++++++++++++++++++++++++++++++++++++

:Data model: `~romancal.datamodels.WfiImage`
:File suffix: _cal

The input to the ``MosaicPipeline`` is a group of calibrated exposures,
e.g. "r0008308002010007027_06311_0019_WFI01_cal.asdf", which contains the
calibrated data for the the exposures. The most convenient way to pass the list of
exposures to be processed with the mosaic level pipeline is to use an association.
Instructions on how to create an input association an be found at :ref:`asn-from-list`.


Outputs
----------

2D Image (MosaicModel)
++++++++++++++++++++++

The resampled data can be found in

:Data model: `~romancal.datamodels.WfiMosaic`
:File suffix: _i2d

Catalog file (MosaicSourceCatalog)
+++++++++++++++++++++++++++++++++++

The catalog data is in

:Data model: `~romancal.datamodels.MosaicSourceCatalog`
:File suffix: _cat

Segmentation Map (SegmentationMapModel)
++++++++++++++++++++++++++++++++++++++++

The segmentation map is

:Data model: `~romancal.datamodels.MosaicSegmentationMapModel`
:File suffix: _segm


Result of applying all the mosaic level pipeline steps up through the
:ref:`source_catalog <source_catalog_step>` step is to produce data background corrected
and cleaned of outliers and resampled to a distortion free grid along with
the source catalog and segmentation map.
The i2d file is 2D image data, with additional attributes for the mosaicing information. The cat
file is an asdf file with the detected sources and the segmenation map is an asdf file
linking the input images to the detected sources.
