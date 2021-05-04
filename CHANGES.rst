0.2.0 (Unreleased)
==================



stpipe
------

- Create stpipe module which provides Roman-specific Step and Pipeline
  subclasses. [#103, #128]

flatfield
---------

- Clean up and improve flatfield step. [#122]

datamodels
----------

- Make necessary changes to use roman_datamodels that is based on the tag approach [#212]

- Added ``RampModel``, ``GLS_RampFitModel``, ``RampFitOutputModel`` and
  schemas. [#110]

- Added ``DQModel`` and schemas. [#81]


0.1.0 (2020-12-11)
==================

datamodels
----------

- First release of romancal. Includes the core metadata and a ``FlatModel``.

- Update date strings in schemas and tests from strings to astropy objects [#32]

-  Update Flat Schema for DQ Array DType [#55]
