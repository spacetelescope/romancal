LINEARITY and INVERSELINEARITY Reference Files
----------------------------------------------

:REFTYPE: LINEARITY, INVERSELINEARITY
:Data models: `~roman_datamodels.datamodels.LinearityRefModel`,
              `~roman_datamodels.datamodels.InverselinearityRefModel`

The LINEARITY reference file contains pixel-by-pixel polynomial coefficients
that map observed DN to linearized DN. The INVERSELINEARITY reference file
has the same format but contains coefficients for the inverse transformation
(linearized DN to observed DN), used internally to simulate individual reads
from observed resultants.

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
detector          model.meta.instrument.detector     WFI
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
