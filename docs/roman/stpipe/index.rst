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
on image and mosaic products as ``meta.cal_logs``. The intended product policy
(see `issue 2510 <https://github.com/spacetelescope/romancal/issues/2510>`_) is:

* **L2 (ELP)**: ELP logs for that exposure. For multi-detector associations,
  each member keeps the log slice from processing that member (plus tweakreg
  messages when alignment completes for that member).
* **L3 (roman_mos)**: mosaic pipeline / resample logs in ``meta.cal_logs``.
  Input L2 logs are retained separately as ``individual_image_cal_logs``.
* **Prompt / forced catalogs**: source-catalog step logs; parent image logs
  may also be included. Upstream logs are not required on segmentation maps.
* **Multiband catalogs**: must not embed L3 production logs. Per-filter
  ``meta.image_metas`` deliberately omits ``cal_logs`` so parquet metadata
  stays free of bulky upstream log streams.

Helpers live in ``romancal.stpipe.cal_logs`` (``set_product_cal_logs``,
``strip_inherited_cal_logs``, ``extract_image_meta``, ``slice_log_records``).

.. automodapi:: romancal.stpipe
