.. _pipelines:

Pipeline Stages
===============

End-to-end calibration of Roman data is divided into levels of
processing:

- The exposure level processing consists of detector-level corrections that are performed on a
  group-by-group basis, followed by ramp fitting. The output of the exposure level
  processing is a count rate image per exposure, or per integration for
  some modes that is aligned to the Gaia reference system.
  Details of this pipeline can be found at:

  - :ref:`exposure_pipeline`

  - Exposure level  processing consists of additional instrument-level and
    observing-mode corrections and calibrations to produce fully calibrated
    exposures. The details differ for imaging and spectroscopic exposures,
    and there are some corrections that are unique to certain instruments or modes.
    Details are at: :ref:`exposure_pipeline`

- The High Level level processing uses overlapping exposures to match the sky background,
  detect aberrant data values and resample the image to produce a single undistorted product.
  Details are at:  :ref:`highlevel_pipeline`

  The table below represents the same information as described above, but
  alphabetically ordered by pipeline class.

+--------------------------------------------+------------------+------------------+
| Pipeline Class                             | Alias            | Used For         |
+============================================+==================+==================+
| `~romancal.pipeline.ExposurePipeline`      | roman_elp        | Exposure Level:  |
+--------------------------------------------+------------------+------------------+
| `~romancal.pipeline.HighLevelPipeline`     | roman_hlp        | High Level:      |
+--------------------------------------------+------------------+------------------+


Pipelines and Exposure Type
===========================

The data from different observing modes needs to be processed with
the proper pipeline levels listed above. The proper pipeline
selection is usually based solely on the exposure type (exposure.type attribute).
