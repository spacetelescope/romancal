Data Product Types
------------------
The following tables contain lists of all data product types.



Exposure Level Data Products
++++++++++++++++++++++++++++

+----------------------------------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+
| Pipeline                                                       | Input                  |  Output(s)               | Base | Units                 | Description                           |
+================================================================+========================+==========================+======+=======================+=======================================+
| :ref:`romancal.pipeline.ExposurePipeline <exposure_pipeline>`  | uncal                  | cal                      | Exp  | DN/s                  | 2-D calibrated data                   |
+----------------------------------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+

Exposure Step Data Products
+++++++++++++++++++++++++++

The following table contain lists of all data product types for exposure pipeline, as given by their file name suffix. The input uncal file and the final cal file
are the only files that are produced in the standard processing. All the other are optional (opt) files that can be produced when
the user is running the pipeline. The input for each optional step is the output of the preceeding step.

+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| Pipeline Step                                  | Input           |  Output(s)               | Model            | Units               | Description                           |
+================================================+=================+==========================+==================+=====================+=======================================+
|                                                |                 | uncal                    | ScienceRawModel  | DN                  | 3-D uncalibrated exposure data        |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`dq_init <dq_init_step>`                  | uncal           | dq_init (opt)            | RampModel        | electron            | 3-D data quality flags applied        |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`saturation <saturation_step>`            |                 | saturation (opt)         | RampModel        | electron            | 3-D data saturated values flagged     |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`reference pixels <refpix>`               |                 | ref_pix (opt)            | RampModel        | electron            | 3-D ref pix corrected data            |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`linearity <linearity_step>`              |                 | linearity (opt)          | RampModel        | electron            | 3-D linearity corrected data          |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`dark_current <dark_current_step>`        |                 | dark_current (opt)       | RampModel        | electron            | 3-D dark current subtracted data      |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`ramp_fitting <ramp_fitting_step>`        |                 | ramp_fit (opt)           | ImageModel       | electron/s          | 2-D slope corrected data              |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`assign_wcs <assign_wcs_step>`            |                 | assign_wcs (opt)         | ImageModel       | electron/s          | 2-D data with gwcs                    |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`flat_field <flatfield_step>`             |                 | flat_field (opt)         | ImageModel       | electron/s          | 2-D QE corrected data                 |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`photom <photom_step>`                    |                 | photom (opt)             | ImageModel       | electron/s          | Add phometric data to header          |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`source_detection <source_detection_step>`|                 | photom (opt)             | ImageModel       | electron/s          | Sources identified in the data        |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`tweakreg <tweakreg_step>`                |                 | tweak_reg (opt)          | ImageModel       | electron/s          | WCS aligned with GAIA                 |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
|                                                |                 | cal                      | ImageModel       | electrons/s         | 2-D calibrated exposure data          |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
