0.19.0 (2025-05-14)
===================

General
-------

- Replace WebbPSF package with STPSF (`#1635
  <https://github.com/spacetelescope/romancal/issues/1635>`_)
- test with latest supported Python version (`#1640
  <https://github.com/spacetelescope/romancal/issues/1640>`_)
- Add resource tracking fixtures for regression tests. (`#1664
  <https://github.com/spacetelescope/romancal/issues/1664>`_)
- Remove meta.group_id and meta.exptype assignments in ModelLibrary. (`#1670
  <https://github.com/spacetelescope/romancal/issues/1670>`_)
- Updates for L1/L2 metadata. (`#1698
  <https://github.com/spacetelescope/romancal/issues/1698>`_)
- Replace root logger usage with module specific loggers. (`#1723
  <https://github.com/spacetelescope/romancal/issues/1723>`_)
- Bump the minimum required version of asdf to 4.1.0 and asdf-astropy to 0.6.0.
  (`#1749 <https://github.com/spacetelescope/romancal/issues/1749>`_)


Documentation
-------------

- Update HLP products for B17 (`#1652
  <https://github.com/spacetelescope/romancal/issues/1652>`_)


Associations
------------

- Adds the programs mk_patchlist and mk_skycellasn to generate level 3
  associations based on skycells in two steps. (`#1657
  <https://github.com/spacetelescope/romancal/issues/1657>`_)
- Update level 3 asn names & add naming utility (`#1753
  <https://github.com/spacetelescope/romancal/issues/1753>`_)


Scripts
-------

- Remove verify_install_requires script which was only for testing and not to
  be packaged. (`#1726
  <https://github.com/spacetelescope/romancal/issues/1726>`_)


``exposure_pipeline``
---------------------

- Write the WfiWcs models after ELP run (`#1680
  <https://github.com/spacetelescope/romancal/issues/1680>`_)


``mosaic_pipeline``
-------------------

- Add ``resample_on_skycell`` option to spec. (`#1642
  <https://github.com/spacetelescope/romancal/issues/1642>`_)


``skycell``
-----------

- match skycells in skymap reference file (`patch_match` -> `skycell.match`)
  (`#1694 <https://github.com/spacetelescope/romancal/issues/1694>`_)


``ramp_fitting`` (WFI-Image, WFI-Prism, WFI-Grism)
--------------------------------------------------

- Remove maker_utils usage. (`#1678
  <https://github.com/spacetelescope/romancal/issues/1678>`_)


``assign_wcs`` (WFI-Image, WFI-Prism, WFI-Grism)
------------------------------------------------

- Remove half pixel offset of s_region. (`#1646
  <https://github.com/spacetelescope/romancal/issues/1646>`_)


``flux``
--------

- Reduce memory usage by computing inplace. (`#1661
  <https://github.com/spacetelescope/romancal/issues/1661>`_)


``tweakreg`` (WFI-Image)
------------------------

- Update create_astrometric_catalog to support ModelLibrary. (`#1646
  <https://github.com/spacetelescope/romancal/issues/1646>`_)


``skymatch``
------------

- Update romancal to use the skymatch code from stcal, instead of having its
  own copy. (`#1465 <https://github.com/spacetelescope/romancal/issues/1465>`_)
- Remove maker_util usage. (`#1678
  <https://github.com/spacetelescope/romancal/issues/1678>`_)


``outlier_detection``
---------------------

- Use skycell wcs from association (if available). Add ``resample_on_skycell``
  option to spec. (`#1642
  <https://github.com/spacetelescope/romancal/issues/1642>`_)
- Change default fillval to nan. (`#1644
  <https://github.com/spacetelescope/romancal/issues/1644>`_)
- Save intermediate products as asdf files with data and wcs keys. (`#1715
  <https://github.com/spacetelescope/romancal/issues/1715>`_)


``resample``
------------

- Use resample code from stcal (`#1634
  <https://github.com/spacetelescope/romancal/issues/1634>`_)
- Use skycell wcs from association (if available). Add ``resample_on_skycell``
  option to spec. (`#1642
  <https://github.com/spacetelescope/romancal/issues/1642>`_)
- Change default fillval to nan. (`#1644
  <https://github.com/spacetelescope/romancal/issues/1644>`_)
- Compute combined wcs from input wcs footprints. (`#1646
  <https://github.com/spacetelescope/romancal/issues/1646>`_)
- Improve runtime by using ``add_model_hook`` from stcal resample. (`#1663
  <https://github.com/spacetelescope/romancal/issues/1663>`_)
- Allow resample to run on files with ``None`` and missing values in metadata.
  (`#1688 <https://github.com/spacetelescope/romancal/issues/1688>`_)


``source_catalog``
------------------

- Changed the flux units in the source catalog from uJy to nJy. (`#1671
  <https://github.com/spacetelescope/romancal/issues/1671>`_)
- Remove use of maker_utils within SourceCatalogStep. (`#1675
  <https://github.com/spacetelescope/romancal/issues/1675>`_)
- Updated column names in source catalog. (`#1686
  <https://github.com/spacetelescope/romancal/issues/1686>`_)
- The source catalogs are now saved in parquet format. (`#1690
  <https://github.com/spacetelescope/romancal/issues/1690>`_)
- Added many new columns to the source catalog. (`#1700
  <https://github.com/spacetelescope/romancal/issues/1700>`_)
- Local background is no longer subtracted from aperture fluxes. (`#1701
  <https://github.com/spacetelescope/romancal/issues/1701>`_)
- Updated column units and added placeholder columns. (`#1731
  <https://github.com/spacetelescope/romancal/issues/1731>`_)
- Update the source catalog documentation. (`#1738
  <https://github.com/spacetelescope/romancal/issues/1738>`_)
- Fix meta.filename in source catalog parquet files. (`#1757
  <https://github.com/spacetelescope/romancal/issues/1757>`_)


``multiband_catalog``
---------------------

- Remove maker_utils usage. (`#1678
  <https://github.com/spacetelescope/romancal/issues/1678>`_)


0.18.0 (2025-02-14)
===================

General
-------

- Remove units from reference file datamodels. (`#1474
  <https://github.com/spacetelescope/romancal/issues/1474>`_)
- remove ``okify_regtests`` script (move to ``ci_watson``) (`#1513
  <https://github.com/spacetelescope/romancal/issues/1513>`_)
- Perform bounding box assignment inline with the ordering that GWCS prefers.
  (`#1527 <https://github.com/spacetelescope/romancal/issues/1527>`_)
- Update romancal to use proper APE 14 API for GWCS interactions. (`#1528
  <https://github.com/spacetelescope/romancal/issues/1528>`_)
- Remove the jump code for the deprecated jump detection for even ramps and
  update the documentation (`#1534
  <https://github.com/spacetelescope/romancal/issues/1534>`_)
- Bump min Python version to 3.11 per SPEC 0. (`#1543
  <https://github.com/spacetelescope/romancal/issues/1543>`_)


Documentation
-------------

- Mention possible need to provide package name to strun when using aliases.
  (`#1476 <https://github.com/spacetelescope/romancal/issues/1476>`_)
- Remove assignment validation example from docs. (`#1504
  <https://github.com/spacetelescope/romancal/issues/1504>`_)


``stpipe``
----------

- Remove Step.__call__ usage (which will be deprecated in stpipe). (`#1499
  <https://github.com/spacetelescope/romancal/issues/1499>`_)


Associations
------------

- Switch association scripts from using ``Main`` class to ``_cli`` function to
  fix return code. (`#1538
  <https://github.com/spacetelescope/romancal/issues/1538>`_)
- This adds additional info to the asn header keyword skycell_wcs_info and
  updates the mosaic pipeline to use
  that information to construct the skycell data from the input exposures.
  (`#1583 <https://github.com/spacetelescope/romancal/issues/1583>`_)
- Fix bug where skycell_wcs_info was double json encoded (`#1592
  <https://github.com/spacetelescope/romancal/issues/1592>`_)


Scripts
-------

- Remove install of missing scripts "schema_editor" and "schemadoc". (`#1538
  <https://github.com/spacetelescope/romancal/issues/1538>`_)
- allow `MosaicModel` in `roman_static_preview` (`#1613
  <https://github.com/spacetelescope/romancal/issues/1613>`_)


``exposure_pipeline``
---------------------

- Fix exposure pipeline handling of all saturated inputs. (`#1525
  <https://github.com/spacetelescope/romancal/issues/1525>`_)
- Update exposure pipeline to use ModelLibrary. (`#1525
  <https://github.com/spacetelescope/romancal/issues/1525>`_)
- Fix description of arguments in docs and add description of fully saturated
  input processing. (`#1593
  <https://github.com/spacetelescope/romancal/issues/1593>`_)


``mosaic_pipeline``
-------------------

- Change the default suffix for mosaic products from _i2d to _coadd (`#1542
  <https://github.com/spacetelescope/romancal/issues/1542>`_)
- Convert calls to ``exit(0)`` within the pipeline into exceptions. (`#1545
  <https://github.com/spacetelescope/romancal/issues/1545>`_)
- Roundtrip L3 wcsinfo especially when skycell specifications are used (`#1585
  <https://github.com/spacetelescope/romancal/issues/1585>`_)


``dq_init`` (WFI-Image, WFI-Prism, WFI-Grism)
---------------------------------------------

- Invoke converter from_tvac_raw to enable processing of TVAC/FPS data (`#1596
  <https://github.com/spacetelescope/romancal/issues/1596>`_)


``saturation`` (WFI-Image, WFI-Prism, WFI-Grism)
------------------------------------------------

- Loosen saturation unit test to allow DO_NOT_USE (`#1571
  <https://github.com/spacetelescope/romancal/issues/1571>`_)
- Add saturation step docs to package index. (`#1593
  <https://github.com/spacetelescope/romancal/issues/1593>`_)


``assign_wcs`` (WFI-Image, WFI-Prism, WFI-Grism)
------------------------------------------------

- Apply velocity aberration correction to the WFI WCS (`#1354
  <https://github.com/spacetelescope/romancal/issues/1354>`_)
- `#1612 <https://github.com/spacetelescope/romancal/issues/1612>`_


``source_detection`` (WFI-Image)
--------------------------------

- Remove SourceDetectionStep (use SourceCatalogStep). (`#1533
  <https://github.com/spacetelescope/romancal/issues/1533>`_)


``tweakreg`` (WFI-Image)
------------------------

- Use PSF astrometry in tweakreg & fix regression test. (`#1578
  <https://github.com/spacetelescope/romancal/issues/1578>`_)


``resample``
------------

- Remove unused arguments from step specification. (`#1593
  <https://github.com/spacetelescope/romancal/issues/1593>`_)


``source_catalog``
------------------

- Use meta.source_catalog and meta.cal_step.source_catalog in source_catalog
  step. (`#1533 <https://github.com/spacetelescope/romancal/issues/1533>`_)
- Create a common mask array for both segmentation and PSF photometry. (`#1574
  <https://github.com/spacetelescope/romancal/issues/1574>`_)
- Add new forced source catalog mode. (`#1611
  <https://github.com/spacetelescope/romancal/issues/1611>`_)


0.23.1 (2025-02-14)
===================

General
-------

- Remove units from reference file datamodels. (`#1474
  <https://github.com/spacetelescope/romancal/issues/1474>`_)
- remove ``okify_regtests`` script (move to ``ci_watson``) (`#1513
  <https://github.com/spacetelescope/romancal/issues/1513>`_)
- Perform bounding box assignment inline with the ordering that GWCS prefers.
  (`#1527 <https://github.com/spacetelescope/romancal/issues/1527>`_)
- Update romancal to use proper APE 14 API for GWCS interactions. (`#1528
  <https://github.com/spacetelescope/romancal/issues/1528>`_)
- Remove the jump code for the deprecated jump detection for even ramps and
  update the documentation (`#1534
  <https://github.com/spacetelescope/romancal/issues/1534>`_)
- Bump min Python version to 3.11 per SPEC 0. (`#1543
  <https://github.com/spacetelescope/romancal/issues/1543>`_)


Documentation
-------------

- Mention possible need to provide package name to strun when using aliases.
  (`#1476 <https://github.com/spacetelescope/romancal/issues/1476>`_)
- Remove assignment validation example from docs. (`#1504
  <https://github.com/spacetelescope/romancal/issues/1504>`_)


``stpipe``
----------

- Remove Step.__call__ usage (which will be deprecated in stpipe). (`#1499
  <https://github.com/spacetelescope/romancal/issues/1499>`_)


Associations
------------

- Switch association scripts from using ``Main`` class to ``_cli`` function to
  fix return code. (`#1538
  <https://github.com/spacetelescope/romancal/issues/1538>`_)
- This adds additional info to the asn header keyword skycell_wcs_info and
  updates the mosaic pipeline to use
  that information to construct the skycell data from the input exposures.
  (`#1583 <https://github.com/spacetelescope/romancal/issues/1583>`_)
- Fix bug where skycell_wcs_info was double json encoded (`#1592
  <https://github.com/spacetelescope/romancal/issues/1592>`_)


Scripts
-------

- Remove install of missing scripts "schema_editor" and "schemadoc". (`#1538
  <https://github.com/spacetelescope/romancal/issues/1538>`_)
- allow `MosaicModel` in `roman_static_preview` (`#1613
  <https://github.com/spacetelescope/romancal/issues/1613>`_)


``exposure_pipeline``
---------------------

- Fix exposure pipeline handling of all saturated inputs. (`#1525
  <https://github.com/spacetelescope/romancal/issues/1525>`_)
- Update exposure pipeline to use ModelLibrary. (`#1525
  <https://github.com/spacetelescope/romancal/issues/1525>`_)
- Fix description of arguments in docs and add description of fully saturated
  input processing. (`#1593
  <https://github.com/spacetelescope/romancal/issues/1593>`_)


``mosaic_pipeline``
-------------------

- Change the default suffix for mosaic products from _i2d to _coadd (`#1542
  <https://github.com/spacetelescope/romancal/issues/1542>`_)
- Convert calls to ``exit(0)`` within the pipeline into exceptions. (`#1545
  <https://github.com/spacetelescope/romancal/issues/1545>`_)
- Roundtrip L3 wcsinfo especially when skycell specifications are used (`#1585
  <https://github.com/spacetelescope/romancal/issues/1585>`_)


``dq_init`` (WFI-Image, WFI-Prism, WFI-Grism)
---------------------------------------------

- Invoke converter from_tvac_raw to enable processing of TVAC/FPS data (`#1596
  <https://github.com/spacetelescope/romancal/issues/1596>`_)


``saturation`` (WFI-Image, WFI-Prism, WFI-Grism)
------------------------------------------------

- Loosen saturation unit test to allow DO_NOT_USE (`#1571
  <https://github.com/spacetelescope/romancal/issues/1571>`_)
- Add saturation step docs to package index. (`#1593
  <https://github.com/spacetelescope/romancal/issues/1593>`_)


``assign_wcs`` (WFI-Image, WFI-Prism, WFI-Grism)
------------------------------------------------

- Apply velocity aberration correction to the WFI WCS (`#1354
  <https://github.com/spacetelescope/romancal/issues/1354>`_)
- `#1612 <https://github.com/spacetelescope/romancal/issues/1612>`_


``source_detection`` (WFI-Image)
--------------------------------

- Remove SourceDetectionStep (use SourceCatalogStep). (`#1533
  <https://github.com/spacetelescope/romancal/issues/1533>`_)


``tweakreg`` (WFI-Image)
------------------------

- Use PSF astrometry in tweakreg & fix regression test. (`#1578
  <https://github.com/spacetelescope/romancal/issues/1578>`_)


``resample``
------------

- Remove unused arguments from step specification. (`#1593
  <https://github.com/spacetelescope/romancal/issues/1593>`_)


``source_catalog``
------------------

- Use meta.source_catalog and meta.cal_step.source_catalog in source_catalog
  step. (`#1533 <https://github.com/spacetelescope/romancal/issues/1533>`_)
- Create a common mask array for both segmentation and PSF photometry. (`#1574
  <https://github.com/spacetelescope/romancal/issues/1574>`_)
- Add new forced source catalog mode. (`#1611
  <https://github.com/spacetelescope/romancal/issues/1611>`_)


0.17.0 (2024-11-15)
===================

General
-------

- Update source catalog file with the tweaked coordinates. (`#1373
  <https://github.com/spacetelescope/romancal/issues/1373>`_)
- move DMS requirement <-> test correlations from ``@metrics_logger()``
  decorators to ``romancal/tests/dms_requirement_tests.json`` (`#1399
  <https://github.com/spacetelescope/romancal/issues/1399>`_)
- Break up long regression tests to avoid needing to okify results twice.
  (`#1426 <https://github.com/spacetelescope/romancal/issues/1426>`_)
- Removed now unused lib.dms. (`#1433
  <https://github.com/spacetelescope/romancal/issues/1433>`_)
- Remove units from romancal. (`#1445
  <https://github.com/spacetelescope/romancal/issues/1445>`_)
- Have pytest clean up some files when it finishes running tests. (`#1446
  <https://github.com/spacetelescope/romancal/issues/1446>`_)
- Fix remaining numpy 2 issues and unpin numpy to allow numpy 2 usage. (`#1447
  <https://github.com/spacetelescope/romancal/issues/1447>`_)
- Give regtest okify results unique subdirectories. (`#1456
  <https://github.com/spacetelescope/romancal/issues/1456>`_)
- Updates to support L1/L2 schema changes. (`#1473
  <https://github.com/spacetelescope/romancal/issues/1473>`_)
- Use stcal to compute s_region keyword. (`#1493
  <https://github.com/spacetelescope/romancal/issues/1493>`_)


Documentation
-------------

- handle changelog entries with ``towncrier`` (`#1375
  <https://github.com/spacetelescope/romancal/issues/1375>`_)
- Update docs to not include default fake values. (`#1419
  <https://github.com/spacetelescope/romancal/issues/1419>`_)


``stpipe``
----------

- Add class_alias for all steps. (`#1509
  <https://github.com/spacetelescope/romancal/issues/1509>`_)


Associations
------------

- add target to asn_from_list command (`#1411
  <https://github.com/spacetelescope/romancal/issues/1411>`_)
- Added code to take an input list of calibrated WFI exposures and creates
  associations based on the
  skycells that they overlap (`#1437
  <https://github.com/spacetelescope/romancal/issues/1437>`_)
- Update skycell_asn docs and add skycell_asn as a script at install time
  (`#1471 <https://github.com/spacetelescope/romancal/issues/1471>`_)
- Updates to the file naming for the products and inputs and adds the
  orientation to the wcs keywords in the asn header (`#1505
  <https://github.com/spacetelescope/romancal/issues/1505>`_)


``mosaic_pipeline``
-------------------

- Allow asn product name to be the output product (`#1394
  <https://github.com/spacetelescope/romancal/issues/1394>`_)


``ramp_fitting`` (WFI-Image, WFI-Prism, WFI-Grism)
--------------------------------------------------

- Drop support for ``ols`` ramp fitting. (`#1398
  <https://github.com/spacetelescope/romancal/issues/1398>`_)


``source_detection`` (WFI-Image)
--------------------------------

- Don't restart loggers during create_gridded_psf_model. (`#1503
  <https://github.com/spacetelescope/romancal/issues/1503>`_)


``tweakreg`` (WFI-Image)
------------------------

- Group by obs_id (`#1448
  <https://github.com/spacetelescope/romancal/issues/1448>`_)
- Updates s_region after running TweakRegStep successfully. (`#1484
  <https://github.com/spacetelescope/romancal/issues/1484>`_)


``outlier_detection``
---------------------

- Update input handling to raise an exception on an invalid input instead of
  issuing a warning and skipping the step. (`#1357
  <https://github.com/spacetelescope/romancal/issues/1357>`_)
- Remove unused arguments to outlier detection. (`#1357
  <https://github.com/spacetelescope/romancal/issues/1357>`_)
- Use stcal common code in outlier detection. (`#1357
  <https://github.com/spacetelescope/romancal/issues/1357>`_)
- Fix bug where on_disk=True could fail due to Quantities not implementing
  tofile. (`#1436 <https://github.com/spacetelescope/romancal/issues/1436>`_)
- Group by obs_id (`#1448
  <https://github.com/spacetelescope/romancal/issues/1448>`_)


``resample``
------------

- Fixed an incompatibility with ``numpy 2.0`` in
  ``resample.resample_utils.build_mask()``. Switched code in
  ``build_driz_weight()`` to use ``astropy`` equivalent of ``build_mask()``.
  Deprecated ``resample.resample_utils.build_mask()``. (`#1383
  <https://github.com/spacetelescope/romancal/issues/1383>`_)
- Group by obs_id (`#1448
  <https://github.com/spacetelescope/romancal/issues/1448>`_)
- Update resample to populate location_name attribute and tests to check for it
  (`#1498 <https://github.com/spacetelescope/romancal/issues/1498>`_)


``source_catalog``
------------------

- The data and err array of the input datamodel to the source_catalog step
  are now copied so that they are left completely unchanged by the step.
  (`#1457 <https://github.com/spacetelescope/romancal/issues/1457>`_)
- Restored flux units in source catalog table. (`#1512
  <https://github.com/spacetelescope/romancal/issues/1512>`_)


``multiband_catalog``
---------------------

- Added a pipeline step to create a multiband catalog from L3 images. (`#1485
  <https://github.com/spacetelescope/romancal/issues/1485>`_)


0.16.3 (2024-08-29)
===================

mosaic_pipeline
---------------

- Only load patch table when needed. [#1367]

source_catalog
--------------

- Populate segmentation image metadata. [#1391]

resample
--------

- Use association product name for output meta.filename by default [#1391]

0.16.2 (2024-08-23)
===================

pipeline
--------

- Added ``suffix`` to the spec of ExposurePipeline with a
  default value of ``cal``. Removed explicit setting of ``suffix``
  so that it can be passed as an argument to ``strun``. [#1378]

0.16.1 (2024-08-13)
===================

- update ``stpipe`` to use ``ModelLibrary`` [#1364]
- update ``stcal`` to use outlier detection [#1364]

0.16.0 (2024-08-13)
===================

Documentation
-------------

- Update RTD to include mosaic data (i2d) description [#1262]

general
-------
- Add regression test for DMS400 and additional tests for ``SkyMatchStep``. [#1358]

- Add regression test for DMS373, mosaic pipeline [#1348]

- Update the exposure pipeline to accept a roman datamodel as input [#1296]

- Update okify script to use GA directory structure [#1282]

- pin numpy to <2 [#1275]

- refactor exposure level pipeline to use asn's and ModelContainer [#1271]

- Add catalog source step to the mosaic pipeline [#1266]

- Rename highlevelpipeline to mosaic pipeline [#1249]

- Replace ``SourceDetectionStep`` with ``SourceCatalogStep`` in ELP. [#1276]

- replace usages of ``copy_arrays`` with ``memmap`` [#1316]

- Replace ModelContainer with ModelLibrary [#1241]

- Updated sky background usage in code and tests to use maker utilities. [#1351]

- Refactor general step input handling to avoid early closing of
  input files to allow using more lazy loading [#1342]



source_catalog
--------------
- Add PSF photometry capability. [#1243]

dq_init
-------
-  Refactor DQInitStep to use the RampModel method of creating ramps. [#1258]

outlier_detection
-----------------

- Set ``single=True`` to use ``many_to_many`` when creating median image. [#1260]

stpipe
------

- Add ``ModelContainer`` support to ``Step._datamodels_open`` to allow
  loading ``pars-*`` files from CRDS. [#1270]


tweakreg
--------
- Integration with ``SourceCatalogStep``: allow usage of results from ``SourceCatalogStep``. [#1276]

resample
--------

- Fix incorrect number of starting planes for context image. [#1355]

mosaic_pipeline
---------------

- Fix construction of skycell WCS.  [#1297]

tweakreg
--------
- Remove unnecessary global variable ALIGN_TO_ABS_REFCAT. [#1314]

- Update default absolute separation for tweakreg.  [#1352]

skymatch
--------
- Populate valid metadata even when then are no overlapping images to
  match [#1360]


0.15.1 (2024-05-15)
===================

- updated `rad` and `roman_datamodels` to `0.20.0`

0.15.0 (2024-05-08)
===================

skymatch
--------
- Update step to always return a ``ModelContainer``. [#1208]

- Fix bug that prevented ``meta.background.subtracted`` from being set with the proper datatype. [#1233]

patch_match
-----------

- Code to determine which patches overlap a given image. [#1161]
- Plotting utility to show image spatial relationship to matched patches and
  candidate patches. [#1204]

tweakreg
--------

- Allow single open Roman datamodels to be used as input to be consistent with expected behavior in ELP. [#1089]

- Update tweakreg regression tests to test astrometric
  performance. Use "clip_accum" for better robustness.  [#1185]

general
-------

- Initial resample to a skycell in the hlp [#1214]

- Add preview files to HLP tests [#1199]

- Allow ``ModelContainer`` to work properly with context manager. [#1147]

- Update the ``dqflags`` to use the ones stored in
  ``roman_datamodels`` [#1099]
- Add script for creating regtest files; consolidate files used for
  some tests. [#1084]

- Update the high level pipeline to use updates in Outlier_detection and tweakreg [#1143]

documentation
-------------

- Fixed datamodels documentation to use correct API. [#1112]

- Improve PSF fitting configuration, background subtraction, grid
  point selection. [#1125]

dq_init
-------

- Copy reference pixels during ``dq_init`` to avoid larger files in later
  processing steps [#1121]

- Allow ``dq_init`` to pass through keys not defined in ``RampModel``
  schema [#1151]

flux
----

- Set flux step status for each input. [#1160]

stpipe
------

- Update ``meta.calibration_software_version`` for results of ``Step`` runs to
  record the version of romancal used to produce the result. [#1194]

- Update ``stpipe.core.finalize_results`` to record the CRDS information
  only if a step uses reference files. [#1201]

- Populate logs for L3 files in addition to L2 files [#1207]

resample
--------

- Update location of ``basic`` attributes. [#1131]

- Allow user to provide DQ flags to use/exclude when creating resampling mask. [#1166]

- Updated Level 3 ``cal_step`` attribute creation. [#1165]

- Fix bug that prevented properly update of the resampled output weight and context arrays. [#1181]

- Update Level 3 output ``basic`` attribute. [#1188]

- Populate the Level 3 wcsinfo [#1182]

- Make rotation matrix 2d for schema validation [#1205]

- Include logs of individual L2 products [#1207]

- Resample members should use actual file names from association file [#1209]

- Populate the l3 product individual_image_meta block [#1216]

outlier_detection
-----------------

- Allow `ModelContainer` as input. [#1092]

- Update location of ``basic`` attributes. [#1131]

- Set ``single=False`` in the call to resample to properly create a median image. [#1146]

ramp_fitting
------------

- Changed image units from e/s to DN/s (and added support for MJy/sr). Added gain reduction to convert to these units. [#1128]

flux
----

- Create FluxStep to apply the flux correction to Level 2 data. [#1120]

source_detection
----------------

- Make PSF fitting the default. [#1185]

source_catalog
--------------

- Added Source Catalog Step. [#1102]

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

- Account for Poisson noise from dark current when fitting ramps. [#1088]

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
