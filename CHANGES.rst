0.3.1 (2021-03-18)
==================

datamodels
----------

- Added sorting to test parameters to preserve order for tests done by paralel pytest workers. [#136]


0.3.0 (2021-03-14)
==================

- Update setup for more strict PEP8 checking [#176]

datamodels
----------

- Updated model subclass code - changed from returning a generator to a set for use with more complicated model selections. [#169]

- Corrected time format in tests to astropy time objects. [#169]

- Cleaned up old tests to better reflect present models. [#169]

- Added check for core metadata inclusion in nonreference files. [#169]


0.2.0 (2021-02-26)
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
