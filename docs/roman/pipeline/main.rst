.. _pipelines:

Pipeline Stages
===============

The data from different observing modes needs to be processed with
the proper pipeline levels listed below. The pipeline, and step selection in a pipeline,
are usually based solely on the exposure type (``exposure.type`` attribute).
End-to-end calibration of Roman data is divided into levels of
processing:

- The Exposure Level Processing (ELP) consists of detector-level corrections applied to
  each resultant, followed by ramp fitting. The output of the exposure level
  processing is a count rate image per exposure, that is aligned to the Gaia reference system.
  The details differ for imaging and spectroscopic exposures and can be found at :ref:`exposure_pipeline`.


- The High Level Processing (HLP) uses overlapping exposures to match the sky background,
  detect aberrant data values and resample the image to produce a single undistorted product.
  Details are at:  :ref:`highlevel_pipeline`

  The table below represents the same information as described above, but
  alphabetically ordered by pipeline class.

+--------------------------------------------+------------------+------------------+
| Pipeline Class                             | Alias            | Used For         |
+============================================+==================+==================+
| `~romancal.pipeline.ExposurePipeline`      | roman_elp        | Exposure Level   |
+--------------------------------------------+------------------+------------------+
| `~romancal.pipeline.HighLevelPipeline`     | roman_hlp        | High Level       |
+--------------------------------------------+------------------+------------------+
