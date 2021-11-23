.. _exposure_pipeline:


exposure_pipeline: Exposure Level Detector Processing
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


============================================== ========= ========= =========
 Step                                          WFI-Image WFI-Prism WFI-Grism
============================================== ========= ========= =========
 :ref:`dq_init <dq_init_step>`                  x          x        x
 :ref:`saturation <saturation_step>`            x          x        x
  `linearity`                                   x          x        x
 :ref:`dark_current <dark_current_step>`        x          x        x
 :ref:`jump <jump_step>`                        x          x        x
 :ref:`ramp_fitting <ramp_fitting_step>`        x          x        x
 :ref:`assign_wcs <assign_wcs_step>`            x          x        x
 :ref:`flatfield <flatfield_step>`              x          x        x
============================================== ========= ========= =========


Arguments
---------
The ``evel2`` pipeline has one optional argument::

  --save_calibrated_ramp  boolean  default=False

If set to ``True``, the pipeline will save intermediate data to a file as it
exists at the end of the :ref:`jump <jump_step>` step (just before ramp fitting). The data
at this stage of the pipeline are still in the form of the original 4D ramps
(ncols x nrows x ngroups x nints) and have had all of the detector-level
correction steps applied to it, including the detection and flagging of
Cosmic-Ray (CR) hits within each ramp (integration). If created, the name of the
intermediate file will be constructed from the root name of the input file, with
the new product type suffix "_ramp" appended,
e.g. "r0008308002010007027_06311_0019_WFI01_ramp.fits".

Inputs
--------

3D raw data
+++++++++++

:Data model: `~romancal.datamodels.RampModel`
:File suffix: _uncal

The input to the ``ExposurePipeline`` is a single raw exposure,
e.g. "r0008308002010007027_06311_0019_WFI01_uncal.fits", which contains the
original raw data from all of the detector readouts in the exposure
(ncols x nrows x ngroups).

Note that in the operational environment, the
input will be in the form of a `~romancal.datamodels.RawScienceModel`, which only
contains the 3D array of detector pixel values, along with some optional
extensions. When such a file is loaded into the pipeline, it is immediately
converted into a `~romancal.datamodels.RampModel`, and has all additional data arrays
for errors and Data Quality flags created and initialized to zero.

Outputs
----------

2D corrected ramp
+++++++++++++++++

:Data model: `~romancal.datamodels.RampModel`
:File suffix: _ramp

Result of applying all pipeline steps up through the :ref:`jump <jump_step>` step,
to produce corrected and CR-flagged 2D ramp data, which will have one less data dimensions
as the input raw 3D data (ncols x nrows x ngroups ).
