LINEARITY Reference File
-------------------------

:REFTYPE: LINEARITY
:Data model: `~roman_datamodels.datamodels.LinearityModel`

The LINEARITY reference file contains pixel-by-pixel polynomial correction
coefficients.

.. include:: ../references_general/linearity_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Attributes for LINEARITY
++++++++++++++++++++++++++++++++++++++
In addition to the standard reference file attributes listed above,
the following attributes are *required* in LINEARITY reference files,
because they are used as CRDS selectors
(see :ref:`linearity_selectors`):

===============  =================================  =============
Attribute         Fully qualified path               Instruments
===============  =================================  =============
detector          model.meta.intstrument.detector    WFI
===============  =================================  =============

Reference File Format
+++++++++++++++++++++
LINEARITY reference files are ASDF format, with 2 data arrays.
The format and content of the file is as follows:

=======  ============ =========================  ===========
Data      Array Type   Dimensions                 Data type
=======  ============ =========================  ===========
coeffs    NDArray      ncols x nrows x ncoeffs    float32
dq        NDArray      ncols x nrows              uint32
=======  ============ =========================  ===========

Each plane of the COEFFS data cube contains the pixel-by-pixel coefficients for
the associated order of the polynomial. There can be any number of planes to
accommodate a polynomial of any order.
