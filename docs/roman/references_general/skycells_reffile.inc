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

Each skycell records the pixel coordinates of its projection tangent
point in ``x_tangent`` and ``y_tangent``, but not the direction of the
pixel x axis, which is needed to interpret them. Files delivered up to
and including ``roman_wfi_skycells_0001.asdf`` place the tangent point
in a mirror-image frame, in which right ascension increases with x. A
later delivery is expected to mirror ``x_tangent`` within each skycell,
giving skycells the conventional orientation, with right ascension
increasing to the left.

`~romancal.skycell.skymap.SkyMap.vparity` infers which convention the
loaded file follows, so both are handled without a code change. Skycell
mosaics built under the two conventions are mirror images of one
another, with the same skycell names and sky footprints.
