.. _area_reffile:

AREA Reference File
-------------------

:REFTYPE: AREA
:Data models: `~roman_datamodels.datamodels.PixelareaRefModel`

The AREA reference file contains pixel area information for a given
detector.

.. include:: ../references_general/area_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for AREA
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in AREA reference files,
because they are used as CRDS selectors
(see :ref:`area_selectors`):

===============   =======================================  ==============
Attribute           Fully qualified path                    Instruments
===============   =======================================  ==============
detector           model.meta.instrument.detector           WFI
optical_element    model.meta.instrument.optical_element    WFI
===============   =======================================  ==============


Reference File Format
+++++++++++++++++++++
AREA reference files are ASDF format, with a 2D data array.
The format and content of the file is as follows
(see `~roman_datamodels.datamodels.PixelareaRefModel`):

=======  ============ ==========================  =============
Data      Array Type   Dimensions                  Data type
=======  ============ ==========================  =============
area      NDArray      4088 x 4088                 float32
=======  ============ ==========================  =============

The area data array contains a 2-D pixel-by-pixel map of relative
pixel areas, normalized to a value of 1.0. The absolute value of the nominal
pixel area is given in the photometry keywords pixelarea_steradians and pixelarea_arcsecsq, in units of steradians and square arcseconds, respectively. The steradian area keyword should have the same values as thepixelareasr keyword in the photometry section of the corresponding PHOTOM
reference file.
