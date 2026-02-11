INTEGRALNONLINEARITY Reference File
-----------------------------------

:REFTYPE: INTEGRALNONLINEARITY
:Data model: `~roman_datamodels.datamodels.IntegralnonlinearityRefModel`

The INTEGRALNONLINEARITY reference file contains per-channel lookup tables
for correcting integral nonlinearity in the analog-to-digital converter.

.. include:: ../references_general/integralnonlinearity_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Attributes for INTEGRALNONLINEARITY
+++++++++++++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file attributes listed above,
the following attributes are *required* in INTEGRALNONLINEARITY reference
files, because they are used as CRDS selectors
(see :ref:`integralnonlinearity_selectors`):

===============  =================================  =============
Attribute         Fully qualified path               Instruments
===============  =================================  =============
detector          model.meta.instrument.detector     WFI
===============  =================================  =============

Reference File Format
+++++++++++++++++++++
INTEGRALNONLINEARITY reference files are ASDF format. The file contains a
lookup value array and correction tables for each of the 32 readout channels.

========================================  ============  ==============  ===========
Data                                      Array Type    Dimensions      Data type
========================================  ============  ==============  ===========
value                                     NDArray       65536           float64
inl_table.science_channel_NN.correction   NDArray       65536           float64
========================================  ============  ==============  ===========

The ``value`` array contains the DN values (0-65535) used as lookup indices.
Each ``science_channel_NN`` (where NN ranges from 01 to 32) contains a
``correction`` array with the INL correction offset for each DN value in that
channel.
