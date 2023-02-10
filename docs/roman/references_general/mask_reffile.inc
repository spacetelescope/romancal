.. _mask_reffile:

MASK Reference File
-------------------

:reftype: MASK
:Data model: `~roman_datamodels.datamodels.MaskRefModel`

The MASK reference file contains pixel-by-pixel DQ flag values that indicate
problem conditions.

.. include:: ../references_general/mask_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for MASK
+++++++++++++++++++++++++++++++
In addition to the standard reference file keyword attributes listed above,
the following keyword attributes are *required* in MASK reference files,
because they are used as CRDS selectors
(see :ref:`mask_selectors`):

===============   ======================================  ==============
Attribute          Fully qualified path                    Instruments
===============   ======================================  ==============
detector           model.meta.instrument.detector          WFI
===============   ======================================  ==============


Reference File Format
+++++++++++++++++++++
MASK reference files are ASDF format, with one data object.
The format and content of the file is as follows:

=======  ============ ==============  =============
Data      Object Type  Dimensions      Data type
=======  ============ ==============  =============
dq        NDArray      4096 x 4096      uint32
=======  ============ ==============  =============

The values in the ``dq`` array give the per-pixel flag conditions that are
to be propagated into the science exposure's ``pixeldq`` array.
The dimensions of the ``dq`` array should be equal to the number of columns
and rows in a full-frame readout of a given detector, including reference
pixels.

The ASDF file contains a single dq array.
