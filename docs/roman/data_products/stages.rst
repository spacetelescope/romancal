Processing Levels and Product Stages
====================================
Here we describe the structure and content of the most frequently used forms of files for
Roman science data products, which are in `ASDF <https://asdf-standard.readthedocs.io/>`_ format. Each type of ASDF
file is the result of serialization of a corresponding :ref:`DataModel<datamodels>`.

Within the various STScI internal data processing and archiving systems that are used for routine processing of
Roman data, there are some different uses of terminology to refer to different levels of processing.
The WFI data is converted into ASDF files by Science Data Formatting (SDF), level 0. SDF also inserts data
from the engineering database and from the proposal database to create the level 1 files. SDF produces one ASDF
file per detector and exposure and these level 1 files are used as input to the Exposure Level Processing. The
output of the exposure level processing is a level 2 file.

+-------------------------------------+-------------------------------------+------------------------------------+
| Data Processing Levels              | User Data Product Stages            | Input data level                   |
+=====================================+=====================================+====================================+
|  Level 0, Raw files                 | Exposure data in an ASDF file       | Uncalibrated                       |
+-------------------------------------+-------------------------------------+------------------------------------+
|  Level 1, Exposure Level Processing | Fully populated ASDF file           | Uncalibrated                       |
+-------------------------------------+-------------------------------------+------------------------------------+
|  Level 2, High Level Processing     | Fully calibrated ASDF file          | calibrated exposures               |
+-------------------------------------+-------------------------------------+------------------------------------+
