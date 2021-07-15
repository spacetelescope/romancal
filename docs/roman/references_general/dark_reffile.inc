.. _dark_reffile:

DARK Reference File
-------------------

:REFTYPE: DARK
:Data models: `~romancal.datamodels.DarkModel`

The DARK reference file contains pixel-by-pixel and frame-by-frame
dark current values for a given detector readout mode.

.. include:: ../references_general/dark_selection.inc

.. include:: ../includes/standard_keywords.inc

.. include:: ../references_general/dark_specific.inc
             
Reference File Format
+++++++++++++++++++++
DARK reference files are ASDF format, with 3 data arrays and one table array.
The format and content of the file is as follows 
(see `~romancal.datamodels.DarkModel`):

=======  ============ ==========================  =============
Data      Array Type   Dimensions                  Data type
=======  ============ ==========================  =============
data      NDArray      4096 x 4096 x ngroups       float32
err       NDArray      4096 x 4096 x ngroups       float32
dq        NDArray      4096 x 4096                 uint16
dq_def    NDArray      2 dims of variable length   various
=======  ============ ==========================  =============

The ASDF file contains a single set of data, err, dq, and dq_def arrays.

#Commented until DQ Docs
#.. include:: ../includes/dq_def.inc
