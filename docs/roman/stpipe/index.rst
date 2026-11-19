======
Stpipe
======

:doc:`Stpipe <stpipe:index>` is used for constructing and running
data pipelines and steps.

.. _stpipe-user-steps:

For Users
=========

See :ref:`running-the-pipeline` for an overview of using the romancal
pipeline.

.. _stpipe-devel-steps:

For Developers
==============

Developers should consult the :doc:`stpipe documentation <stpipe:index>`
for details about constructing new steps, pipelines and using logging.

.. automodapi:: romancal.stpipe

Calibration logs
================

``RomanStep.finalize_result`` stores the producing step or pipeline log
stream on ``meta.cal_logs`` for image, mosaic, source-catalog, and
segmentation products. Catalog and segmentation products share the same
catalog-step records. Mosaic products do not retain per-input
``individual_image_cal_logs``; mosaic production messages live in
``meta.cal_logs``.
