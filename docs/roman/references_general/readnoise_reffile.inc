READNOISE reference file
------------------------

:REFTYPE: READNOISE
:Data model: `~roman_datamodels.datamodels.ReadnoiseRefModel`

The READNOISE reference file contains a pixel-by-pixel map of readnoise, which is
used in estimating the expected noise in each pixel.

.. include:: ../references_general/readnoise_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keyword Attributes for READNOISE
++++++++++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file keyword attributes listed above,
the following keyword attributess are *required* in READNOISE reference files,
because they are used as CRDS selectors
(see :ref:`readnoise_selectors`):

===============  ============================  ===========
Keyword Attr.    Fully qualified path          Instruments
===============  ============================  ===========
detector         model.meta.wfi_mode.detector  WFI
===============  ============================  ===========

Reference File Format
+++++++++++++++++++++
READNOISE reference files are asdf files with a single data array.

=======  ============ ==============  =============
Data      Array Type   Dimensions       Data type
=======  ============ ==============  =============
data      NDArray      4096 x 4096      float32
=======  ============ ==============  =============

The ASDF file contains a single data array.
