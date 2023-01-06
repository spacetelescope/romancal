.. _saturation_reffile:

SATURATION Reference File
-------------------------

:REFTYPE: SATURATION
:Data model: `~roman_datamodels.datamodels.SaturationRefModel`

The SATURATION reference file contains pixel-by-pixel saturation threshold
values.

.. include:: ../references_general/saturation_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for SATURATION
+++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in SATURATION reference files,
because they are used as CRDS selectors
(see :ref:`saturation_selectors`):

=========  ==============================  ==========
Keyword    Data Model Name                 Instrument
=========  ==============================  ==========
detector   model.meta.instrument.detector  WFI
=========  ==============================  ==========

Reference File Format
+++++++++++++++++++++
SATURATION reference files are ASDF format, with two data objects.
The format and content of the file is as follows:

=======  ============ ==============  =============
Data      Object Type  Dimensions      Data type
=======  ============ ==============  =============
data      NDArray      4096 x 4096      float32
dq        NDArray      4096 x 4096      uint32
=======  ============ ==============  =============

The values in the ``data`` array give the saturation threshold in units of DN
for each pixel.

The ASDF file contains two data arrays.
