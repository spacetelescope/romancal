.. _photom_reffile:

PHOTOM Reference File
---------------------

:REFTYPE: PHOTOM
:Data models: `~roman_datamodels.datamodels.WfiImgPhotomRefModel`

The PHOTOM reference file contains conversion factors for putting
pixel values into physical units.

.. include:: ../references_general/photom_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for PHOTOM
+++++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in PHOTOM reference files. The first (detector) is
needed because it is used as a CRDS selector. The second (optical_element) is used to
select the appropriate set of photometric parameters.
(see :ref:`photom_selectors`):

===============   =======================================  ==============
Attribute          Fully qualified path                     Instruments
===============   =======================================  ==============
detector           model.meta.instrument.detector           WFI
optical_element    model.meta.instrument.optical_element    WFI
===============   =======================================  ==============


Tabular PHOTOM Reference File Format
++++++++++++++++++++++++++++++++++++
PHOTOM reference files are ASDF format, with data in the phot_table attribute.
The format and content of the file is as follows
(see `~romancal.datamodels.WfiImgPhotomRefModel`):

Data is stored in a 2D table, with optical elements for the row names:

+------------+---------------------------------------------------------------------+
| Instrument | Row names                                                           |
+============+=====================================================================+
| WFI        | F062, F087, F106, F129, F146, F158, F184, F213, GRISM, PRISM, DARK  |
+------------+---------------------------------------------------------------------+


And the variable attributes for the columns (with data type):

+------------+----------------+-----------+------------+------------------------+
| Instrument | Column name    | Data type | Dimensions | Units                  |
+============+================+===========+============+========================+
| WFI        | photmjsr       | quantity  | scalar     | MJy/steradian          |
+            +----------------+-----------+------------+------------------------+
|            | uncertainty    | quantity  | scalar     | MJy/steradian          |
+            +----------------+-----------+------------+------------------------+
|            | pixelareasr    | quantity  | scalar     | steradian              |
+------------+----------------+-----------+------------+------------------------+

The pixelareasr variable attribute gives the average pixel area in units of steradians.
