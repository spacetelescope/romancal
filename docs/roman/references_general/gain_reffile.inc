GAIN reference file
-------------------

:REFTYPE: GAIN
:Data model: `~roman_datamodels.datamodels.GainRefModel`

The GAIN reference file contains a pixel-by-pixel gain map, which can be
used to convert pixel values from units of DN to electrons. The gain
values are assumed to be in units of e/DN.

.. include:: ../references_general/gain_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keyword Attributes for GAIN
+++++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keyword attributes listed above,
the following keyword attributess are *required* in GAIN reference files,
because they are used as CRDS selectors
(see :ref:`gain_selectors`):

===============  ============================  ===========
Keyword Attr.    Fully qualified path          Instruments
===============  ============================  ===========
detector         model.meta.wfi_mode.detector  WFI
===============  ============================  ===========

Reference File Format
+++++++++++++++++++++
GAIN reference files are asdf files with a single data array.

=======  ============ ==============  =============
Data      Array Type   Dimensions       Data type
=======  ============ ==============  =============
data      NDArray      4096 x 4096      float32
=======  ============ ==============  =============

The ASDF file contains a single data arrays.
