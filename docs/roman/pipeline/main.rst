.. _pipelines:

Pipeline Stages
===============

End-to-end calibration of Roman data is divided into stages of
processing:

- Stage 2 consists of detector-level corrections that are performed on a
  group-by-group basis, followed by ramp fitting. The output of stage 1
  processing is a count rate image per exposure, or per integration for
  some modes. Details of this pipeline can be found at:

  - :ref:`exposure_pipeline`

  - Stage 2 processing consists of additional instrument-level and
    observing-mode corrections and calibrations to produce fully calibrated
    exposures. The details differ for imaging and spectroscopic exposures,
    and there are some corrections that are unique to certain instruments or modes.
    Details are at: TBD



  The table below represents the same information as described above, but
  alphabetically ordered by pipeline class.

+------------------------------------------+-----------------+-----------+
| Pipeline Class                           | Alias           | Used For  |
+=========================================+==================+===========+
| `~romancal.pipeline.ExposurePipeline`    | roman_elp       | Stage 2:  |
+-----------------------------------------+------------------+-----------+


Pipelines and Exposure Type
===========================

The data from different observing modes needs to be processed with
the proper pipeline stages listed above. The proper pipeline
selection is usually based solely on the exposure type (exposure.type attribute).
