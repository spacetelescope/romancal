.. _flat_reffile:

FLAT Reference File
-------------------

:REFTYPE: FLAT
:Data model: `~romancal.datamodels.FlatModel`

The FLAT reference file contains pixel-by-pixel detector response values.
It is used.

.. include:: ../references_general/flat_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for FLAT
+++++++++++++++++++++++++++++++
In addition to the standard reference file attributes listed above,
the following attributes are *required* in FLAT reference files,
because they are used as CRDS selectors
(see :ref:`flat_selectors`):

===============   ======================================  ==============
Attribute           Fully qualified path                   Instruments
===============   ======================================  ==============
detector           model.meta.instrument.detector          WFI
exptype            model.meta.exposure.type
optical_element    model.meta.instrument.optical_element
===============   ======================================  ==============

Reference File Format
+++++++++++++++++++++

FLAT reference files are ASDF format, with 3 data arrays.
The format and content of the file is as follows:

=======  ============ ==============  =============
Data      Array Type   Dimensions       Data type
=======  ============ ==============  =============
data      NDArray      4096 x 4096      float32
err       NDArray      4096 x 4096      float32
dq        NDArray      4096 x 4096      uint16
=======  ============ ==============  =============

The ASDF file contains a single set of data, err, and dq arrays.
