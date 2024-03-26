Processing Levels and Product Stages
====================================
Here we describe the structure and content of the most frequently used forms of files for
Roman science data products, which are in `ASDF <https://asdf-standard.readthedocs.io/>`_ format. Each type of ASDF
file is the result of serialization of a corresponding :ref:`DataModel<datamodels>`. All data models are
defined by their `schemas <https://rad.readthedocs.io/latest/>`_.

Within the various STScI internal data processing and archiving systems that are used for routine processing of
Roman data, there are some different uses of terminology to refer to different levels of processing.
The WFI data is converted into ASDF files by Science Data Formatting (SDF), level 0. SDF also inserts data
from the engineering database and from the proposal database to create the level 1 files. SDF produces one ASDF
file per detector and exposure and these level 1 files are used as input to the Exposure Level Processing. The
output of the exposure level processing is a level 2 file.

+-------------------------------------+---------------------------------------------+----------------------------+
| Data Processing Levels              | User Data Product Stages                    | MAST Product               |
+=====================================+=============================================+============================+
|  Level 0, Science telemetry         | Not available to users                      | Not available to users     |
+-------------------------------------+---------------------------------------------+----------------------------+
|  Level 1, Uncalibrated files        | Level 1 - ASDF file, fully populated by SDF | 1 - Uncalibrated exposures |
+-------------------------------------+---------------------------------------------+----------------------------+
|  Exposure Level Processing          | Level 2 - Fully calibrated exposures        | 2 - Calibrated exposures   |
+-------------------------------------+---------------------------------------------+----------------------------+
|  High Level Processing              | Level 3 - Combined data                     | 3 - Combined data          |
+-------------------------------------+---------------------------------------------+----------------------------+
|  High Level Processing              | Level 4 - Catalogs, segmentation images     | 4 - Combined data          |
+-------------------------------------+---------------------------------------------+----------------------------+


Level 1 and Level 2 products packaged the data from individual detectors in a single ASDF file. Level 3
are data resampled on a sky cell. The definitions of the sky cells is TBD. The definition of Level 3 and
Level 4 products is being finalized.
