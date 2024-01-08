.. _outlier_design:

OutlierDetectionStep
--------------------

This module provides the sole interface to all methods of performing outlier detection
on Roman observations. The outlier detection algorithm used for WFI data is implemented
in :py:class:`~romancal.outlier_detection.outlier_detection.OutlierDetection`
and described in :ref:`outlier-detection-imaging`.

.. note::
    Whether the data are being provided in an `association file`_ or as a list of ASDF filenames,
    they must always be wrapped with a
    :py:class:`~romancal.datamodels.container.ModelContainer`, which will handle and
    read in the input properly.

.. _association file: https://jwst-pipeline.readthedocs.io/en/latest/jwst/associations/asn_from_list.html

On successful completion, this step will return the input models with DQ arrays updated
with flags for identified outliers.


.. automodapi:: romancal.outlier_detection.outlier_detection_step
