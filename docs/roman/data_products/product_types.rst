Data Product Types
------------------

The following tables contain lists of all data product types.



Exposure Level Data Products
++++++++++++++++++++++++++++

+------------------------------------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+
| Pipeline                                                         | Input                  |  Output(s)               | Base | Units                 | Description                           |
+==================================================================+========================+==========================+======+=======================+=======================================+
| :ref:`romancal.pipeline.ExposurePipeline <exposure_pipeline>`    | uncal                  | cal                      | Exp  | DN/s                  | 2-D calibrated data                   |
+------------------------------------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+
| :ref:`romancal.pipeline.HighLevelPipeline <highlevel_pipeline>`  | cal                    | i2d                      | Exp  | DN/s                  | 2-D calibrated data                   |
+------------------------------------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+


Exposure Pipeline Steps And Data Products
+++++++++++++++++++++++++++++++++++++++++

The following table contains lists of all data product types for the Exposure pipeline, as given by their file name suffix. The input uncal file and the final cal file
are the only files that are produced in the standard processing. All the other are optional (opt) files that can be produced when
the user is running the pipeline. The input for each optional step is the output of the preceding step.

+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| Pipeline Step                                  | Input           |  Output suffix           | Data Model       | Units               | Description                           |
+================================================+=================+==========================+==================+=====================+=======================================+
|                                                |                 | uncal                    | ScienceRawModel  | DN                  | 3-D uncalibrated exposure data        |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`dq_init <dq_init_step>`                  | uncal           | dqinit (opt)             | RampModel        | DN                  | 3-D data quality flags applied        |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`saturation <saturation_step>`            |                 | saturation (opt)         | RampModel        | DN                  | 3-D data saturated values flagged     |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`refpix <refpix>`                         |                 | refpix (opt)             | RampModel        | DN                  | 3-D ref pix corrected data            |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`linearity <linearity_step>`              |                 | linearity (opt)          | RampModel        | DN                  | 3-D linearity corrected data          |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`dark_current <dark_current_step>`        |                 | darkcurrent (opt)        | RampModel        | DN                  | 3-D dark current subtracted data      |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`ramp_fitting <ramp_fitting_step>`        |                 | rampfit (opt)            | ImageModel       | DN/s                | 2-D slope corrected data              |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`assign_wcs <assign_wcs_step>`            |                 | assignwcs (opt)          | ImageModel       | DN/s                | 2-D data with gwcs                    |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`flatfield <flatfield_step>`              |                 | flat (opt)               | ImageModel       | DN/s                | 2-D QE corrected data                 |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`photom <photom_step>`                    |                 | photom (opt)             | ImageModel       | DN/s                | Add phometric data to header          |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`source_detection <source_detection_step>`|                 | sourcedetectionstep (opt)| ImageModel       | DN/s                | Sources identified in the data        |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
| :ref:`tweakreg <tweakreg_step>`                |                 | tweakregstep (opt)       | ImageModel       | DN/s                | WCS aligned with GAIA                 |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+
|                                                |                 | cal                      | ImageModel       | DN/s                | 2-D calibrated exposure data          |
+------------------------------------------------+-----------------+--------------------------+------------------+---------------------+---------------------------------------+



High Level Processing Steps And Data Products
+++++++++++++++++++++++++++++++++++++++++++++

The following table contain lists of all data product types for the HighLevel Processing (HLP) Pipeline, as given by their file name suffix.
The input to the HLP is an association file (in JSON format), the output is a combined image.
All the other are optional (opt) files that can be produced when
the user is running the pipeline. The input for each optional step is the output of the preceding step.

+---------------------------------------------------+-----------------+------------------------------+------------------+---------------------+---------------------------------------+
| Pipeline Step                                     | Input           |  Output suffix               | Data Model       | Units               | Description                           |
+===================================================+=================+==============================+==================+=====================+=======================================+
|                                                   |                 | asn                          |                  |                     |                                       |
+---------------------------------------------------+-----------------+------------------------------+------------------+---------------------+---------------------------------------+
| :ref:`sky_match <skymatch_step>`                  | asn             | skymatch (opt)               | ModelContainer   | MJy/sr              | A list of _cal files                  |
+---------------------------------------------------+-----------------+------------------------------+------------------+---------------------+---------------------------------------+
| :ref:`outlier_detection <outlier_detection_step>` |                 | outlier_detection_step (opt) | ModelContainer   | MJy/sr              | A list of _cal files                  |
+---------------------------------------------------+-----------------+------------------------------+------------------+---------------------+---------------------------------------+
| :ref:`resample <resample_step>`                   |                 | resamplestep (opt)           | ModelContainer   | MJy/sr              | A list of _cal files                  |
+---------------------------------------------------+-----------------+------------------------------+------------------+---------------------+---------------------------------------+
|                                                   |                 | i2d                          | MosaicModel      | MJy/sr              | A 2D resampled image                  |
+---------------------------------------------------+-----------------+------------------------------+------------------+---------------------+---------------------------------------+
