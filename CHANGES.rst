0.7.1 (2022-05-19)
==================

general
-------

- Update regression tests with new data, remove skips for flat fielding tests, and code cleanup [#504]

jump
----

- Enable multiprocessing in jump detection step. [#503]

linearity
---------

- Account for possible zero frame in linearity [#506]

saturation
----------

- Updated the saturation step due to an update in STCAL. [#500]

0.7.0 (2022-05-13)
==================

general
-------

- Update regression tests with new data conforming to the latest datamodels [#496]

- Remove copy of arrays that are not needed to help manage memory [#487]

Documentation
-------------

- Add documentation for error propagation in ramp fitting and flat field [#476]

- Add documentation for DNS build 0.5, e.g. reference array trimming [#457]

- Updated documentation for the photom step and removed the area reference
  documentation. [#488]

- Added documentation for Distortion reference files. [#493]


linearity
---------

-  Linearity correction now supports NaN's in the reference file. [#484]

  photom
------

- Photom updated to skip updating photometric converstions for spectral data [#498]

- Added photom correction step and unit tests. [#469]

- Added SOC test for absolute photometric calibration. Tweak logging in photom step. [#479]


0.6.0 (2022-03-02)
==================

general
-------

- Update the regression test for new datamodels and suffixes. [#442]

- Updated PEP 8 checks to be more comprehensive. [#417]

- Added regression tests for linearity correction. [#394]

- Added regression tests for dark_current subtraction. [#392]

- Updated tests to utilize new maker function code. [#395]

- Border reference pixel arrays (and their dq) are copied in ``dq_init``.
  They are trimmed from the science data (and err/dq) in ``ramp_fit``. [#435]

Documentation
-------------

 - Add documentation on using info and search with Roman datamodels [#432]

 - Add the suffixes used in the pipeline if steps.<step>.save_results is set [#415]

 - Update references_general.rst to remove TBD and add DQ flag information. [#396]

 - Initial romancal documentation for using datamodels. [#391]

 - Added documentation for PHOTOM and Area reference files, which required placeholder documentation for the photom step. In addition, I fixed an improper object in dark documentation. [#452]

dark
----

 - Updated dark current step to use stcal. Created tests for the updated step. [#420]

 - Fixed dark subtraction output crash. [#423]


jump
----

 - Update Jump regression test parameters to reduce test time [#411]

 - Update code to suppress output from the jump step if not requested [#399]

Pipeline
________
 - Migrate JWST suffix infrastructure to the Roman Exposure Pipeline [#425]


0.5.0 (2021-12-13)
==================

general
-------

- Added regression tests for SOC-604. [#381]

- Added regression tests for SOC-622. [#385]


linearity
---------

- Implemented linearity correction using stcal. [#360]

assign_wcs
----------

- Added ``assign_wcs`` step to romancal. [#361]

flat
----

- Added check in flat field step to skip spectroscopic observations. Added test. [#366]

jump
----

- Updated filenames in regression test script [#351]

- Updates to add the suffix _flat to the step output [#349]

- Updates for unit tests to use stcal [#322]

- Fix to jump_step to save the update pixel and group dq arrays. [#319]

- Updated code for ``jump`` step using ``stcal``. [#309]

- Added simple regression test. [#315]

- Updated temp readnoise file in jump tests to include required exposure keywords. [#333]

ramp_fitting
------------

- Update ramp_fitting regression test output file names [#369]

- Implemented ramp_fitting using stcal. [#276]

saturation
----------

- Implement saturation correction using stcal, roman_datamodels and romancal.stpipe [#348]

- Updated RTD to include saturation reference files. [#350]

stpipe
------

 - Record step/pipeline logs in ImageModel.cal_logs array. [#352]

0.4.2 (2021-09-13)
==================

general
-------

- Corrected artifactory path from romancal-pipeline to roman-pipeline. [#295]

0.4.1 (2021-09-02)
==================

general
-------

- Updated requirements-sdp.txt for release.


0.4.0 (2021-09-01)
==================

general
-------

- Added regressions tests for ``dq_ini``t utilizing ``mask`` file in CRDS. [#290]

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
==================

general
-------
- Added grism to the CRDS tests [# 225]


0.3.0 (2021-05-28)
==================

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
