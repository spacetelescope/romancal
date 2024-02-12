0.14.1 (unreleased)
==================



0.14.0 (2024-02-12)
===================

general
-------

- Updated the ``compare_asdf`` diff reports to include descriptive information
  about what is being compared. [#1044]

dq_init
-------

- Add the ability to copy resultantdq from a SDF science raw model to the new rampmodel created by dq_init [#1085]

outlier_detection
-----------------

- Add outlier detection step documentation. [#1042]
- Add outlier detection unit tests. [#1058]
- Add additional documentation of the scale and snr parameters. [#1058]
- Updated information for the ``scale`` and ``snr`` parameters in the ``outlier_detection`` step docs. [#1062]

jump detection
--------------

- Added uneven ramp-jump detection docs. [#1035]


documentation
-------------

- Remove ``sphinx-asdf`` requirement, fix issue where menu does not scroll. [#1063]

- Update jump step docs [#1035]

- added user documentation for ``roman_static_preview`` script [#1046]


ramp_fitting
------------

- Add default WCS when constructing image model from ramp model [#1072]

resample
--------

- Update resample step to handle the L3 meta data [#1057]


general
-------

- Update elp steps to check for CRDS not returning a reference file [#1055]

- Fix bug where ``compare_asdf`` failed to detect ``DataModel`` type differences. [#1066]


0.13.0 (2023-11-28)
===================

outlier_detection
-----------------

- Implemented ``outlier-detection step``. [#981]

associations
------------

- Add FOV associations to the  code  [#931]

dark
----

- Removed ``err`` array from dark current tests. [#938]

general
-------

- Update elp pipeline code to capture a list from tweakreg [#985]

- Add code to run the steps needed for the high level processing (roman_hlp) [#980]

- Update pipeline code to correct cal_step and suffixes [#971]

- Update pipeline code to run through tweakreg with single files and associations [#960]

- Update regression tests with new data and update ramp fitting tests to use ols_cas22 [#911]

- Fix bug with ``ModelContainer.get_crds_parameters`` being a property not a method [#846]

- Fix random seed bug in PSF fitting methods [#862]

- Fix regression tests for PSF fitting methods [#872]

- Fix regression test ``compare_asdf`` function replacing use of
  ``asdf.commands.diff`` with ``deepdiff`` and add ``deepdiff`` as
  a test dependency [#868]

- Add ``astropy.table.Table`` support to ``compare_asdf`` [#915]

- Use tolerance for more comparisons in ``compare_asdf`` [#917]

- Use array comparison options (including ``nan`` equality) when
  comparing ``WCS`` objects during ``compare_asdf`` [#941]

- Fix dynamic importing issue with the ``ddtrace`` package. [#1024]

ramp_fitting
------------

- Inititial implementation of the Uneven Ramp fitting [#779]

- Fix opening mode for references to be read-only [#854]

- Make uneven ramp fitting the default [#877]

- Update Ramp fitting code to support the ``stcal`` changes to the ramp fitting
  interface which were necessary to support jump detection on uneven ramps [#933]

- Add uneven ramp fitting documentation [#944]

- Enable jump detection within the Cas22 ramp fitting be default, and add
  regression tests for it. [#991]

- Implement next round of SOC verification tests for uneven ramps [#970]

refpix
------

- Update cal_step, add suffix and add to the exposure pipeline [#890]

- Enable apodized FFT interpolation by default. [#1017]

resample
--------

- Implement resampling step. [#787]

- Use resampled exposure time images to compute image exposure times.  [#959]

scripts
-------

- added ``roman_static_preview`` script to generate static previews of ASDF images [#953]

- fixed ``asn_from_list`` script [#972]

source_detection
----------------

- Support for PSF fitting (optional) for accurate centroids. [#841, #984]

- Save source catalog to a structured array. [#987]

stpipe
------

- Remove checks on CI in production code [#955]

tweakreg
--------

- Fix a bug due to which source catalog may contain sources
  outside of the bounding box. [#947]

0.12.0 (2023-08-18)
===================

source_detection
----------------
- Skip the step if the data is not imaging mode. [#798]

tweakreg
--------
- Skip the step if the data is not imaging mode [#798]

- Add regression test for TweakReg. [#707]

- WCS fit results are now available in meta.wcs_fit_results. [#714]

documentation
-------------
- Update info strings in the pipeline to provide uniform syntax [#721]

- Updated wording about ELP and HLP in the Associations documentation for RTD

- Updated the primary branch referenced in CONTRIBUTING to be main

- Updated reference pixel correction documentation to include discretization bias discussion. [#716]

skymatch
--------
- Added SkyMatchStep to pipeline [#687]

- Registered SkyMatchStep in stpipe. [#770]

jump
----
- Accept and ignore additional return values from stcal detect_jumps [#723]

ramp_fitting
------------
- Update unit tests for stcal 1.4.0 [#725]

- Adjust ramp slopes and associated unceratinties for gain. [#804]

refpix
------

- Add initial reference pixel correction step implementation. [#704]

saturation
----------

- Add read_pattern argument to flag_saturated_pixels. [#836]

general
-------

- Add metrics_logger to the regression tests [#831]

- Update pipeline logic for saturation checks [#824]

- Update the pipeline code to process all the uncal files in an association [#802]

- `ModelContainer` supports slice and dice. [#710]

- Add `ModelContainer` to `romancal.datamodels`. [#710]

- Move ``is_assocation`` from ``roman_datamodels`` to ``romancal``. [#719]

- Update ``romancal`` to use altered API for ``maker_utils``. [#717]

- Require stcal >= 1.4 [#723]

- Fix search for docs. [#768]

- Remove ``aws`` install option. [#767]

- Bump minimum ``asdf`` version to ``2.15.0``. [#777]

- Remove unused extras (``ephem``, ``lint``) from build configuration and regression testing [#784]

- Make all random number generation for tests both seeded and use the same random
  number generation system. [#771]

- Make steps operate in place rather than copying.  [#774]

- Fix devdeps Jenkins job. [#795]

- Remove use of the deprecated ``pkg_resources`` module from ``setuptools``. [#829]

- Add ``dev`` install option. [#835]

- Add PSF photometry methods [#794]

0.11.0 (2023-05-31)
===================

tweakreg
--------

- Added tmpdir to the unit tests for test files [#702]

- Added logic to handle cases where an absolute catalog cannot be created. [#698]

associations
------------

- Initial association code for GBTDS observations [#661]

Documentation
-------------

- Update dq flags to include "GW_AFFECTED_DATA"  flag [#699]

general
-------
- Updated datamodel maker utility imports. [#654]

- Update non-VOunits to using ``astropy.units``. [#658]

- update minimum version of ``asdf`` to ``2.14.2`` and ``jsonschema`` to ``4.0.1`` and added minimum dependency checks to CI [#664]

- Remove use of ``pytest-openfiles`` [#666]

- Remove the ``codecov`` dependency [#677]

- Remove explicit dependence on ``stdatamodels``. [#676]

- Drop support for Python 3.8 [#694]

source_detection
----------------
- Bug fix to ensure that the returned result is a copy of the input datamodel. [#700]

- Added SourceDetection Step to pipeline [#608]

- Added option of fixed random seed for unit tests to avoid intermittent failures from randomness. [#668]

- Fix source detection object instantiation. [#669]

- Small bug fix to ensure that output catalogs are not attached to the file when save_catalogs=False [#684]

outlier_detection
-----------------
- Added an empty outlier detection step to the pipeline, as well as a simple test and documentation. [#689]

astrometric_utils
-----------------
- Added option to provide epoch so that the coordinates are corrected by proper motion. [#686]


0.10.0 (2023-02-21)
===================

general
-------
- Adds explicit test for PSF keywords are present in the  cal files. [#648]

- Add ``pre-commit`` configuration to repository. [#622]

- Use ``isort`` and ``black`` to format code, also upgrade all string
  formats using ``flynt``. [#645]

- Update the suffix for the stored filename to match the filename [#609]

- DQ step flags science data affected by guide window read [#599]

- Fix deprecation warnings introduced by ``pytest`` ``7.2`` ahead of ``8.0`` [#597]

- Implemented support for quantities in reference files. Updated unit tests for these changes. [#624]

associations
------------

- Initial association code with asn_from_list and some basic rules [#642]


jump
----

- Update jump units to roman_datamodels from astropy units [#646]

- Update default input CR thresholds to give reasonable results [#625]

- Added support for Quantities for data arrays. [#616]

tweakreg
--------
- First implementation of TweakRegStep into the pipeline [#643]


0.9.0 (2022-11-14)
==================

general
-------

- New Roman's RTD page layout [#596]

- pin ``numpy`` to ``>=1.20`` [#592]
- replace ``flake8`` with ``ruff`` [#570]


jump
----

- Changes for new keywords (currently unused by Roman) to control snowball and shower flagging in jump detection. [#593]

photom
------

- Updates so that the default suffix is used for spectroscopic data. [#594]

- Change photom step to forcibly set the photometric keywords to ``None`` for spectroscopic data. [#591]

tests
-----

- refactor `tox` environment factors and structure GitHub Actions into dependent workflow [#551]

0.8.1 (2022-08-23)
==================

- pin ``asdf`` above ``2.12.1`` to fix issue with `jsonschema` release [#562]

- pin `roman_datamodels` to newest feature version [#563]

0.8.0 (2022-08-12)
==================

assign_wcs
----------

- Add distortion transform to assign_wcs step. [#510]

Documentation
-------------

- include information about the distortion reference file used in the ``assign_wcs`` step [#542]

flat
----

- Removed try/except condition on Flat Reference file CRDS lookup. [#528]

general
-------

- Update pipeline steps to define the default suffix when saving the step results [#521]
- Simplified reference file name and model storage in dq and flat steps. [#514]

- Update CI workflows to cache test environments and depend upon style and security checks [#511]
- Release ``numpy`` version requirement [#544]
- Moved build configuration from ``setup.cfg`` to ``pyproject.toml`` to support PEP621 [#512]
- Added support for STCAL handing of fully saturated data in both the pipeline and rampfit step. Added a unit test for the rampfit changes and a regression test for the pipeline chages. [#541]

- Update `stpipe` requirement to `>=0.4.2` [#545]

- Fix input_filename when DataModel is input to ExposurePipeline [#553]

- Populate 'ref_file' section in meta after step is run. [#492]

- pin ``asdf`` above ``2.12.1`` to fix issues with unit and regression tests [#562]

photom
------

- Adds explicit test that photometric keywords are preserved for spectroscopic data. [#513]

- Changed optical element W146 to F146. [#552]


ramp_fitting
------------

- Added multiprocessing ramp test. Fixed ols ramp fit. Updated ramp_fit to add photometry to image file generation. [#523]

tests
-----

- Updated tests to account for the change in dimensionality of the err variable in ramp datamodel. [#520]
- Added SOC tests to check for information available in Level 2 images to correct for pixel geometric distortion. [#549]

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

Documentation
-------------

- Add documentation for error propagation in ramp fitting and flat field [#476]

- Add documentation for DNS build 0.5, e.g. reference array trimming [#457]

- Updated documentation for the photom step and removed the area reference
  documentation. [#488]

- Added documentation for Distortion reference files. [#493]

- Updated wording about ELP and HLP in the Associations documentation for RTD

- Updated the primary branch referenced in CONTRIBUTING to be main


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

- Added regressions tests for ``dq_init`` utilizing ``mask`` file in CRDS. [#290]

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
