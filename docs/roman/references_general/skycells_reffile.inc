.. _skycells_reffile:

SkyCells Reference File
-----------------------

:REFTYPE: SKYCELLS
:Data models: `~roman_datamodels.datamodels.SkycellsRefModel`

The SKYCELLS reference file describes the tesselation of
projection regions over the sphere, as well as skycells within
each projection region.

Handedness
++++++++++

The ``meta.vparity`` keyword gives the direction of the pixel x axis.
Files delivered up to and including ``roman_wfi_skycells_0001.asdf``
lack it and use a mirror-image frame (``vparity = +1``), in which right
ascension increases with x.  Subsequent reference files set
``vparity = -1``, giving skycells the conventional orientation, with
right ascension increasing to the left.  A skycell mosaic built under
the new convention is the mirror image, in x, of one built on the same
footprint under the old convention.
