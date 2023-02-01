:orphan:

.. _distortion_reffile:

DISTORTION Reference File
-------------------------

:REFTYPE: DISTORTION
:Data model: `~roman_datamodels.datamodels.DistortionRefModel`

The distortion reference file contains a combination of astropy models,
representing the transform from detector to the telescope V2, V3 system.
The following convention was adopted:

- the input x and y are 0-based coordinates in the DMS system;
- the output is in the V2, V3 system;
- the center of the first pixel is (0, 0), so the first pixel goes from -0.5 to 0.5;
- the origin of the transform is taken to be (0, 0).
  Note that while a different origin can be used for some transforms the relevant
  offset should first be prepended to the distortion transform to account for the change
  in origin of the coordinate frame.

Internally the WCS pipeline works with 0-based coordinates.

.. include:: ../references_general/distortion_selection.inc

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
DISTORTION reference files are ASDF format, and contain an astropy model object.
The format and content of the file is as follows
(see `~roman_datamodels.datamodels.DistortionRefModel`):

==================================  ========================
Data                                 Data type
==================================  ========================
coordinate_distortion_transform      astropy.modeling.Model
==================================  ========================

The ASDF file contains a single astropy model object.

:model: Transform from detector to an intermediate frame (instrument dependent).
