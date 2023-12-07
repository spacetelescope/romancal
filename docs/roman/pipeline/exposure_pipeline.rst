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
 :ref:`linearity <linearity_step>`                  |check|    |check|  |check|
 :ref:`dark_current <dark_current_step>`            |check|    |check|  |check|
 :ref:`jump <jump_step>`                            |check|    |check|  |check|
 :ref:`ramp_fitting <ramp_fitting_step>`            |check|    |check|  |check|
 :ref:`assign_wcs <assign_wcs_step>`                |check|    |check|  |check|
 :ref:`flatfield <flatfield_step>`                  |check|
 :ref:`photom <photom_step>`                        |check|
 :ref:`source_detection <source_detection_step>`   |check|
 :ref:`tweakreg <tweakreg_step>`                    |check|
================================================== ========= ========= =========


Arguments
---------
The ``exposure`` pipeline has an optional argument::

  --use_ramp_jump_detection  boolean  default=True

When set to ``True``, the pipeline will perform :ref:`jump <jump_step>`  detection as a part of the ramp
fitting  step. The data at this stage of the pipeline are still in the form of the original
3D ramps ( ngroups x ncols x nrows ) and have had all of the detector-level
correction steps applied to it, up to but not including the detection and flagging of
Cosmic-Ray (CR) hits within each ramp (integration). For this case the  :ref:`jump <jump_step>`
module in :ref:`ramp_fitting <ramp_fitting_step>` will update the dq array with the CR hits (jumps) that
are identified in the step.


Inputs
--------

3D raw data
+++++++++++

:Data model: `~romancal.datamodels.RampModel`
:File suffix: _uncal

The input to the ``ExposurePipeline`` is a single raw exposure,
e.g. "r0008308002010007027_06311_0019_WFI01_uncal.asdf", which contains the
original raw data from all of the detector readouts in the exposure
( ngroups x ncols x nrows ).

Note that in the operational environment, the
input will be in the form of a `~romancal.datamodels.RawScienceModel`, which only
contains the 3D array of detector pixel values, along with some optional
extensions. When such a file is loaded into the pipeline, it is immediately
converted into a `~romancal.datamodels.RampModel`, and has all additional data arrays
for errors and Data Quality flags created and initialized to zero.

Outputs
----------

2D Image model
++++++++++++++

:Data model: `~romancal.datamodels.ImageModel`
:File suffix: _cal

Result of applying all pipeline steps up through the
:ref:`tweakreg <tweakreg_step>` step is to produce calibrated data with the image wcs
aligned to Gaia, and is 2D image data, which will have one less data dimensions as the input
raw 3D data ( ngroups x ncols x nrows ). In addition to being a 2-dimensional
image the output from the pipeline has the :ref:`reference pixels <refpix>`
removed from the edges of the science array.
