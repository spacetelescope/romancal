.. _exposure_pipeline:


Exposure Level Processing
=====================================================

:Class: `romancal.pipeline.ExposurePipeline`
:Alias: exposure_pipeline

The ``ExposurePipeline`` applies detector-level corrections to given exposure
types (imaging, prism, and grism.). It is applied to one
exposure at a time.
It is sometimes referred to as "ramps-to-slopes" processing, because the input
raw data are in the form of ramps containing accumulating counts from the
non-destructive detector readouts and the output is a corrected countrate
(slope) image.

The list of steps applied by the ``ExposurePipeline`` pipeline is shown in the
table below.

.. |check| unicode:: U+2713 .. checkmark

================================================== ========= ========= =========
 Step                                              WFI-Image WFI-Prism WFI-Grism
================================================== ========= ========= =========
 :ref:`dq_init <dq_init_step>`                      |check|    |check|  |check|
 :ref:`saturation <saturation_step>`                |check|    |check|  |check|
 :ref:`refpix <refpix>`                             |check|    |check|  |check|
 :ref:`linearity <linearity_step>`                  |check|    |check|  |check|
 :ref:`dark_current <dark_current_step>`            |check|    |check|  |check|
 :ref:`ramp_fitting <ramp_fitting_step>`            |check|    |check|  |check|
 :ref:`assign_wcs <assign_wcs_step>`                |check|    |check|  |check|
 :ref:`flatfield <flatfield_step>`                  |check|
 :ref:`photom <photom_step>`                        |check|
 :ref:`source_catalog <source_catalog_step>`        |check|
 :ref:`tweakreg <tweakreg_step>`                    |check|
================================================== ========= ========= =========


Arguments
---------
The ``exposure`` pipeline has no arguments

Inputs
------

3D raw data
+++++++++++

:Data model: `~romancal.datamodels.RampModel`
:File suffix: _uncal

The input to the ``ExposurePipeline`` can be a single raw exposure,
e.g. "r0008308002010007027_0019_wfi01_uncal.asdf", which contains the
original raw data from all of the detector readouts in the exposure
( ngroups x ncols x nrows ). The raw data may also be input using an association file.

If the ``ExposurePipeline`` is given a single file the final alignment to Gaia will be done
with the sources found in the exposure. If multiple exposures exist in the association file
then the final alignment will use all the sources found in the exposures
(see :ref:`tweakreg <tweakreg_step>`).

Note that in the operational environment, the
input will be in the form of a `~romancal.datamodels.RawScienceModel`, which only
contains the 3D array of detector pixel values, along with some optional
extensions. When such a file is loaded into the pipeline, it is immediately
converted into a `~romancal.datamodels.RampModel`, and has all additional data arrays
for errors and Data Quality flags created and initialized.

When the ``ExposurePipeline`` processes a fully saturated input (all pixels flagged as saturated).
The corresponding output image will:

- contain all 0 data arrays
- contain all 0 variance arrays
- not be processed by steps beyond saturation

A single fully saturated input will also cause :ref:`tweakreg <tweakreg_step>` to be skipped
for all input images.

Outputs
-------

2D Image model
++++++++++++++

:Data model: `~romancal.datamodels.ImageModel`
:File suffix: _cal

Catalog file (SourceCatalog)
+++++++++++++++++++++++++++++++++++

The catalog data is in

:Data model: `astropy.table.Table`
:File suffix: _cat

Segmentation Map (SegmentationMapModel)
++++++++++++++++++++++++++++++++++++++++

The segmentation map is

:Data model: `~romancal.datamodels.MosaicSegmentationMapModel`
:File suffix: _segm

Result of applying all pipeline steps up through the
:ref:`tweakreg <tweakreg_step>` step is to produce calibrated data with the image WCS
aligned to Gaia, and is 2D image data, which will have one less data dimensions as the input
raw 3D data. In addition to being a 2-dimensional
image the output from the pipeline has the :ref:`reference pixels <refpix>`
removed from the edges of the science array and saved as additional 3D arrays. The
source catalog and segmentation map from the individual exposues is also saved.

WFI Level 1/Level 2 WCS (WfiWcsModel)
+++++++++++++++++++++++++++++++++++++

:Data model: `~romancal.datamodels.WfiWcsModel`
:File suffix: _wcs

Contains a copy of the final, GAIA-aligned, Level 2 Generalized World Coordinate
System (GWCS) information along with a modified Level 1 GWCS which accounts for
the border pixels. The Level 1 GWCS can be used directly with the related Level
1 product.
