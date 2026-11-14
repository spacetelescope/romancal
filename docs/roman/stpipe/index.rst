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

Calibration logs on products
----------------------------

``RomanStep.finalize_result`` records the current step or pipeline log stream
on image and mosaic products as ``meta.cal_logs``.

* **L2 images**: full ELP log stream for the run that produced the file.
* **L3 mosaics**: mosaic/resample pipeline logs in ``meta.cal_logs``. Input L2
  logs are kept separately on ``individual_image_cal_logs``.
* **Source catalogs (parquet)**: typically do not carry a first-class
  ``cal_logs`` field; multiband catalogs omit L3 ``cal_logs`` from
  ``meta.image_metas`` so they are not flattened into parquet metadata.
* **Segmentation maps**: no ``cal_logs`` field.

.. automodapi:: romancal.stpipe
