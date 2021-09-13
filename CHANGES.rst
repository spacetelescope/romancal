0.4.2 (2021-09-13)


general
----------

- Corrected artifactory path from romancal-pipeline to roman-pipeline. [#295]


0.4.1 (2021-09-02)

general
----------

- updated requirements-sdp.txt for release


0.4.0 (2021-09-01)

general
----------

- Added regressions tests for dq_init utilizing mask file sin CRDS. [#290]

- Updates for requirements & pip changes [#286]

- Added test for crds flat file temporal matching (SOC-636.1). [#283]

- Updates for readthedocs [#260]

- Added DQ support. [#262]

- Added stcal as dependency on romancal [#255]

- Locked romancal library dependency version RDM (0.1.2). [#246]

- Update roman_datamodels, stcal, and stpipe to resolve issues with recent
  pip releases. [#284]

Documentation
-------------

- Updated README weblinks.[#241]

- Added documentation for dark current reference files. [#232]

- Added documentation for gain step. [#231]



0.3.1 (2021-06-02)
=======

general
-------

- Added grism to the CRDS tests [# 225]

0.3.0 (2021-05-28)
=======

datamodels
----------

- Added sorting to test parameters to preserve order for tests done by parallel pytest workers. [#136]

- Update setup.cfg to match JWST warnings & error list and initial pass for code fixes. (#188)

general
-------
- Added grism to the regression tests [# 222]

- Update README and CHANGES.rst [#195]

- Added sorting to test parameters to preserve order for tests done by parallel
  pytest workers. [#136]

- Update setup for more strict PEP8 checking [#176]

- Added documentation for rmask files. [#181]

datamodels
----------

- Make necessary changes to use roman_datamodels that is based on the tag approach [#212]

- Add cal_step added to datamodels [#177]

- Updated model subclass code - changed from returning a generator to a set
  for use with more complicated model selections. [#169]

- Corrected time format in tests to astropy time objects. [#169]

- Cleaned up old tests to better reflect present models. [#169]

- Added check for core metadata inclusion in non-reference files. [#169]

- Add Photom Schema [#200]

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

- Add unit tests for the dark current subtraction step [#168]

- Add dark current subtraction step for use with WFI data [#146]

- Add datamodel and schema for mask files [#143]

- Update output_ext in the base Step class to .asdf from .fits [#127]


=======

- Added ``RampModel``, ``GLS_RampFitModel``, ``RampFitOutputModel`` and
  schemas. [#110]

- Update core schema with latest filter information [#97]

- Add the variable arrays to the schema & datamodel for Image files [#93]

- Add Roman Readnoise model [#90]

- Add Gain Model Schema [#82]

- Added ``DQModel`` and schemas. [#81]


0.1.0 (2020-12-11)
==================

datamodels
----------

- First release of romancal. Includes the core metadata and a ``FlatModel``.

- Update date strings in schemas and tests from strings to astropy objects [#32]

- Add Ramp Model Schema [#56]

- Update Flat Schema for DQ Array DType [#55]

- Add exptype information for roman data [#41]

- Use Astropy Time Objects in date and Useafter [#32]

- Add level 1 schema file for Wide Field Imaging model [#31]

- Create a Data Models sub-package for Roman [#17]

- Use the ASDF pytest plugin to validate the datamodels schemas [#6]
