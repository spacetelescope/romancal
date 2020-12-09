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
In addition to the standard reference file keywords listed above,
the following keywords are *required* in FLAT reference files,
because they are used as CRDS selectors
(see :ref:`flat_selectors`):

=========  ==============================  =============================
Keyword    Data Model Name                 Instruments
=========  ==============================  =============================
DETECTOR   model.meta.instrument.detector  All
EXP_TYPE   model.meta.exposure.type        
FILTER     model.meta.instrument.filter
BAND       model.meta.instrument.band      
READPATT   model.meta.exposure.readpatt    
GRATING    model.meta.instrument.grating   
=========  ==============================  =============================

Reference File Format
+++++++++++++++++++++
FLAT reference files are ASDF format, with 3 IMAGE arrays
and 1 table. The format and content of the file is as follows:

=======  ========  =====  ==============  ==============================
EXTNAME  XTENSION  NAXIS  Dimensions      Data type
=======  ========  =====  ==============  ==============================
SCI      IMAGE       2    ncols x nrows   float
ERR      IMAGE       2    ncols x nrows   float
DQ       IMAGE       2    ncols x nrows   integer
DQ_DEF   TABLE       NA   4 columns       integer, integer, string, string
=======  ========  =====  ==============  ==============================

The ``DQ_DEF`` table extension lists the bit assignments for the flag conditions
used in the DQ array.

.. include:: ../includes/dq_def.inc

The ASDF file contains a single set of SCI, ERR, and DQ arrays, and a DQ_DEF
table.  Image dimensions should be 2048x2048.

