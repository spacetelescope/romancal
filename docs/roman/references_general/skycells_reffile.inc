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
in a mirror-image frame, in which right ascension increases with x.
``roman_wfi_skycells_0002.asdf`` mirrors ``x_tangent`` within each
skycell, giving skycells the conventional orientation, with right
ascension increasing to the left.

`~romancal.skycell.skymap.SkyMap.vparity` infers which convention the
loaded file follows, so both are handled without a code change. A
skycell mosaic built under the new convention is the mirror image, in
x, of one built on the same footprint under the old convention.

``roman_wfi_skycells_0002.asdf`` also renames the skycells, so that the
``x`` and ``y`` indices in a name follow the pixel x and y axes. The
footprint named ``xAAyBB`` in the old file is named ``x(100-BB)yAA`` in
the new one; for example, ``270p65x70y49`` became ``270p65x51y70``. The
set of footprints is unchanged.
