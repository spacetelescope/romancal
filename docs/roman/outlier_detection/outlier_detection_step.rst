.. _outlier_design:

OutlierDetectionStep
-----------------------------------------

This module provides the sole interface to all methods of performing outlier
detection on JWST observations.  The ``outlier_detection`` step supports multiple
algorithms and determines the appropriate algorithm for the type of observation
being processed.


This step uses the following logic to apply the appropriate algorithm to the
input data:

* Interpret inputs (ASN table, ModelContainer or CubeModel)
  to identify all input observations to be processed

* Read in type of exposures in input by interpreting ``meta.exposure.type`` from inputs

* Read in parameters set by user

* Select outlier detection algorithm based on exposure type

  - **Images**: like those taken with NIRCam, will use
    :py:class:`~romancal.outlier_detection.outlier_detection.OutlierDetection` as described
    in :ref:`outlier-detection-imaging`

* Instantiate and run outlier detection class determined for the exposure type
  using parameter values interpreted from inputs.

* Return input models with DQ arrays updated with flags for identified outliers


.. automodapi:: romancal.outlier_detection.outlier_detection_step
